//  Copyright 2007-2010 Portland State University, University of Wisconsin-Madison
//  Author: Robert Scheller, Melissa Lucash

using Landis.Core;
using Landis.SpatialModeling;
using Landis.Library.LeafBiomassCohorts;
using System.Collections.Generic;
using Landis.Library.Climate;

namespace Landis.Extension.Succession.NECN_Hydro
{
    /// <summary>
    /// Utility methods.
    /// </summary>
    public class Century
    {
        public static int Year;
        public static int Month;
        public static int MonthCnt;

        /// <summary>
        /// Grows all cohorts at a site for a specified number of years.
        /// Litter is decomposed following the Century model.
        /// </summary>
        public static ISiteCohorts Run(ActiveSite site,
                                       int         years,
                                       bool        isSuccessionTimeStep)
        {
            
            ISiteCohorts siteCohorts = SiteVars.Cohorts[site];
            IEcoregion ecoregion = PlugIn.ModelCore.Ecoregion[site];

            for (int y = 0; y < years; ++y) {

                Year = y + 1;
                
                if (Climate.Future_MonthlyData.ContainsKey(PlugIn.FutureClimateBaseYear + y + PlugIn.ModelCore.CurrentTime - years))
                    ClimateRegionData.AnnualWeather[ecoregion] = Climate.Future_MonthlyData[PlugIn.FutureClimateBaseYear + y - years + PlugIn.ModelCore.CurrentTime][ecoregion.Index];

                SiteVars.ResetAnnualValues(site);

                if (y == 0 && SiteVars.FireSeverity != null && SiteVars.FireSeverity[site] > 0)
                    FireEffects.ReduceLayers(SiteVars.FireSeverity[site], site);

                // Next, Grow and Decompose each month
                int[] months = new int[12]{6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5};

                if(OtherData.CalibrateMode)
                    //months = new int[12]{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}; This output will not match normal mode due to differences in initialization
                    months = new int[12] { 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5 };

                PlugIn.AnnualWaterBalance = 0;

                for (MonthCnt = 0; MonthCnt < 12; MonthCnt++)
                {
                    // Calculate mineral N fractions based on coarse root biomass.  Only need to do once per year.
                    if (MonthCnt == 0)
                    {
                        AvailableN.CalculateMineralNfraction(site);
                    }

                    Month = months[MonthCnt];

                    SiteVars.MonthlyAGNPPcarbon[site][Month] = 0.0;
                    SiteVars.MonthlyBGNPPcarbon[site][Month] = 0.0;
                    SiteVars.MonthlyNEE[site][Month] = 0.0;
                    SiteVars.MonthlyResp[site][Month] = 0.0;
                    SiteVars.MonthlyStreamN[site][Month] = 0.0;
                    SiteVars.SourceSink[site].Carbon = 0.0;
                    SiteVars.TotalWoodBiomass[site] = Century.ComputeWoodBiomass(site);
                                   
                    double ppt = ClimateRegionData.AnnualWeather[ecoregion].MonthlyPrecip[Century.Month];

                    double monthlyNdeposition;
                    if (PlugIn.AtmosNintercept != -1 && PlugIn.AtmosNslope != -1)
                    {
                        monthlyNdeposition = PlugIn.AtmosNintercept + (PlugIn.AtmosNslope * ppt);
                    }
                    else 
                    {
                        monthlyNdeposition = ClimateRegionData.AnnualWeather[ecoregion].MonthlyNDeposition[Century.Month];
                    }

                    ClimateRegionData.MonthlyNDeposition[ecoregion][Month] = monthlyNdeposition;
                    ClimateRegionData.AnnualNDeposition[ecoregion] += monthlyNdeposition;
                    SiteVars.MineralN[site] += monthlyNdeposition;

                    double liveBiomass = (double) ComputeLivingBiomass(siteCohorts);
                    double baseFlow, stormFlow, AET;
                    SoilWater.Run(y, Month, liveBiomass, site, out baseFlow, out stormFlow, out AET);

                    PlugIn.AnnualWaterBalance += ppt - AET;

                    // Calculate N allocation for each cohort
                    AvailableN.SetMineralNallocation(site);

                    if (MonthCnt==11)
                        siteCohorts.Grow(site, (y == years && isSuccessionTimeStep), true);
                    else
                        siteCohorts.Grow(site, (y == years && isSuccessionTimeStep), false);

                    WoodLayer.Decompose(site);
                    LitterLayer.Decompose(site);
                    SoilLayer.Decompose(site);

                    // Volatilization loss as a function of the mineral N which remains after uptake by plants.  
                    // ML added a correction factor for wetlands since their denitrification rate is double that of wetlands
                    // based on a review paper by Seitziner 2006.

                    double volatilize = (SiteVars.MineralN[site] * PlugIn.DenitrificationRate); 

                    SiteVars.MineralN[site] -= volatilize;
                    SiteVars.SourceSink[site].Nitrogen += volatilize;
                    SiteVars.Nvol[site] += volatilize;

                    SoilWater.Leach(site, baseFlow, stormFlow);

                    SiteVars.MonthlyNEE[site][Month] -= SiteVars.MonthlyAGNPPcarbon[site][Month];
                    SiteVars.MonthlyNEE[site][Month] -= SiteVars.MonthlyBGNPPcarbon[site][Month];
                    SiteVars.MonthlyNEE[site][Month] += SiteVars.SourceSink[site].Carbon;
                    SiteVars.FineFuels[site] = (SiteVars.SurfaceStructural[site].Carbon + SiteVars.SurfaceMetabolic[site].Carbon) * 2.0;

                }
            }

            ComputeTotalCohortCN(site, siteCohorts);

            return siteCohorts;
        }

        //---------------------------------------------------------------------

        public static int ComputeLivingBiomass(ISiteCohorts cohorts)
        {
            int total = 0;
            if (cohorts != null)
                foreach (ISpeciesCohorts speciesCohorts in cohorts)
                    foreach (ICohort cohort in speciesCohorts)
                        total += (int) (cohort.WoodBiomass + cohort.LeafBiomass);
                    //total += ComputeBiomass(speciesCohorts);
            return total;
        }

        //---------------------------------------------------------------------

        public static double ComputeWoodBiomass(ActiveSite site)
        {
            double woodBiomass = 0;
            if (SiteVars.Cohorts[site] != null)
                foreach (ISpeciesCohorts speciesCohorts in SiteVars.Cohorts[site])
                    foreach (ICohort cohort in speciesCohorts)
                        woodBiomass += cohort.WoodBiomass;
            return woodBiomass;
        }

        //---------------------------------------------------------------------
        public static void ComputeTotalCohortCN(ActiveSite site, ISiteCohorts cohorts)
        {
            SiteVars.CohortLeafC[site] = 0;
            SiteVars.CohortFRootC[site] = 0;
            SiteVars.CohortLeafN[site] = 0;
            SiteVars.CohortFRootN[site] = 0;
            SiteVars.CohortWoodC[site] = 0;
            SiteVars.CohortCRootC[site] = 0;
            SiteVars.CohortWoodN[site] = 0;
            SiteVars.CohortCRootN[site] = 0;

            if (cohorts != null) 
                foreach (ISpeciesCohorts speciesCohorts in cohorts)
                    foreach (ICohort cohort in speciesCohorts)
                        CalculateCohortCN(site, cohort);
        }

        /// <summary>
        /// Summarize cohort C&N for output.
        /// </summary>
        private static void CalculateCohortCN(ActiveSite site, ICohort cohort)
        {
            ISpecies species = cohort.Species;

            double leafC = cohort.LeafBiomass * 0.47;
            double woodC = cohort.WoodBiomass * 0.47;

            double fRootC = Roots.CalculateFineRoot(cohort, leafC);
            double cRootC = Roots.CalculateCoarseRoot(cohort, woodC);

            double leafN  = leafC / SpeciesData.LeafCN[species];
            double woodN  = woodC / SpeciesData.WoodCN[species];
            double cRootN = cRootC / SpeciesData.CoarseRootCN[species];
            double fRootN = fRootC / SpeciesData.FineRootCN[species];

            SiteVars.CohortLeafC[site]  += leafC;
            SiteVars.CohortFRootC[site] += fRootC;
            SiteVars.CohortLeafN[site]  += leafN;
            SiteVars.CohortFRootN[site] += fRootN;
            SiteVars.CohortWoodC[site]  += woodC;
            SiteVars.CohortCRootC[site] += cRootC;
            SiteVars.CohortWoodN[site]  += woodN ;
            SiteVars.CohortCRootN[site] += cRootN;
        }

    }
}
