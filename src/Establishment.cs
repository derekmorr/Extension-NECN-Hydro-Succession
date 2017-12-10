//  Copyright 2007-2016 Portland State University
//  Author: Robert Scheller

using Landis.Core;
using Landis.SpatialModeling;
using Edu.Wisc.Forest.Flel.Util;
using Landis.Library.Succession;

using System;
using System.IO;
using Landis.Library.Climate;

namespace Landis.Extension.Succession.NECN_Hydro
{
    public class Establishment
    {

        private static StreamWriter log;
        private static double[,] avgSoilMoisturelimit = new double[PlugIn.ModelCore.Species.Count, PlugIn.ModelCore.Ecoregions.Count]; 
        private static double[,] avgMATlimit = new double[PlugIn.ModelCore.Species.Count, PlugIn.ModelCore.Ecoregions.Count]; 
        private static double[,] avgJanuaryTlimit = new double[PlugIn.ModelCore.Species.Count, PlugIn.ModelCore.Ecoregions.Count]; 
        private static double[,] avgPest = new double[PlugIn.ModelCore.Species.Count, PlugIn.ModelCore.Ecoregions.Count]; 


        public static void InitializeLogFile()
        {
            string logFileName   = "NECN_Hydro-prob-establish-log.csv"; 
            PlugIn.ModelCore.UI.WriteLine("   Opening a NECN_Hydro log file \"{0}\" ...", logFileName);
            try {
                log = Landis.Data.CreateTextFile(logFileName);
            }
            catch (Exception err) {
                string mesg = string.Format("{0}", err.Message);
                throw new System.ApplicationException(mesg);
            }
            
            log.AutoFlush = true;
            log.WriteLine("Time, Species, ClimateRegion, AvgTempMult, AvgMinJanTempMult, AvgSoilMoistureMult, AvgProbEst");
        }

        public static double Calculate(ISpecies species, ActiveSite site)
        {
            IEcoregion climateRegion = PlugIn.ModelCore.Ecoregion[site];

            double tempMultiplier = 0.0;
            double soilMultiplier = 0.0;
            double minJanTempMultiplier = 0.0;
            double establishProbability = 0.0;

            AnnualClimate_Monthly ecoClimate = ClimateRegionData.AnnualWeather[climateRegion];

            if (ecoClimate == null)
                throw new System.ApplicationException("Error in Establishment: CLIMATE NULL.");

            double ecoDryDays = SiteVars.DryDays[site];
            soilMultiplier = SoilMoistureMultiplier(ecoClimate, species, ecoDryDays);
            tempMultiplier = BotkinDegreeDayMultiplier(ecoClimate, species);
            minJanTempMultiplier = MinJanuaryTempModifier(ecoClimate, species);

            // Liebig's Law of the Minimum is applied to the four multipliers for each year:
            double minMultiplier = System.Math.Min(tempMultiplier, soilMultiplier);
            minMultiplier = System.Math.Min(minJanTempMultiplier, minMultiplier);

            establishProbability += minMultiplier;
            establishProbability *= PlugIn.ProbEstablishAdjust;

            avgSoilMoisturelimit[species.Index, climateRegion.Index] += soilMultiplier;
            avgMATlimit[species.Index, climateRegion.Index] += tempMultiplier;
            avgJanuaryTlimit[species.Index, climateRegion.Index] += minJanTempMultiplier;
            avgPest[species.Index, climateRegion.Index] += establishProbability;

            return establishProbability;
        }

        public static void LogEstablishment()
        {
            foreach (ISpecies species in PlugIn.ModelCore.Species)
            {
                foreach (IEcoregion ecoregion in PlugIn.ModelCore.Ecoregions)
                {
                    if (!ecoregion.Active || ClimateRegionData.ActiveSiteCount[ecoregion] < 1)
                        continue;

                        log.Write("{0}, {1}, {2},", PlugIn.ModelCore.CurrentTime, species.Name, ecoregion.Name);
                        log.Write("{0:0.00},", (avgMATlimit[species.Index, ecoregion.Index] / (double)ClimateRegionData.ActiveSiteCount[ecoregion]));
                        log.Write("{0:0.00},", (avgJanuaryTlimit[species.Index, ecoregion.Index] / (double)ClimateRegionData.ActiveSiteCount[ecoregion]));
                        log.Write("{0:0.00},", (avgSoilMoisturelimit[species.Index, ecoregion.Index] / (double)ClimateRegionData.ActiveSiteCount[ecoregion]));
                        log.WriteLine("{0:0.00}", (avgPest[species.Index, ecoregion.Index] / (double)ClimateRegionData.ActiveSiteCount[ecoregion]));
                }
            }

        avgSoilMoisturelimit = new double[PlugIn.ModelCore.Species.Count, PlugIn.ModelCore.Ecoregions.Count];
        avgMATlimit = new double[PlugIn.ModelCore.Species.Count, PlugIn.ModelCore.Ecoregions.Count];
        avgJanuaryTlimit = new double[PlugIn.ModelCore.Species.Count, PlugIn.ModelCore.Ecoregions.Count];
        avgPest = new double[PlugIn.ModelCore.Species.Count, PlugIn.ModelCore.Ecoregions.Count];
    }


    //---------------------------------------------------------------------------
    private static double SoilMoistureMultiplier(AnnualClimate weather, ISpecies species, double dryDays)

        {
            double sppAllowableDrought = SpeciesData.MaxDrought[species];
            double growDays = 0.0;
            double maxDrought;
            double Soil_Moist_GF = 0.0;

            growDays = weather.EndGrowing - weather.BeginGrowing + 1.0;
            if (growDays < 2.0)
            {
                PlugIn.ModelCore.UI.WriteLine("Begin Grow = {0}, End Grow = {1}", weather.BeginGrowing, weather.EndGrowing);
                throw new System.ApplicationException("Error: Too few growing days.");
            }
            //Calc species soil moisture multipliers
            maxDrought = sppAllowableDrought * growDays;
            
            if (maxDrought < dryDays) 
            {
                Soil_Moist_GF = 0.0;
            }
            else
            {
                Soil_Moist_GF = System.Math.Sqrt((maxDrought - dryDays) / maxDrought);
            }
            return Soil_Moist_GF;
        }
        
        //---------------------------------------------------------------------------
        private static double BotkinDegreeDayMultiplier(AnnualClimate weather, ISpecies species)
        {

            //Calc species degree day multipliers  
            //Botkin et al. 1972. J. Ecol. 60:849 - 87
            
            double max_Grow_Deg_Days = SpeciesData.GDDmax[species]; 
            double min_Grow_Deg_Days = SpeciesData.GDDmin[species];
            
            double Deg_Day_GF = 0.0;
            double Deg_Days = (double) weather.GrowingDegreeDays; 
            double totalGDD = max_Grow_Deg_Days - min_Grow_Deg_Days;
            
            Deg_Day_GF = (4.0 * (Deg_Days - min_Grow_Deg_Days) * 
                  (max_Grow_Deg_Days - Deg_Days)) / (totalGDD * totalGDD);
            
           if (Deg_Day_GF < 0) Deg_Day_GF = 0.0;
           
           return Deg_Day_GF;
        }
        
        //---------------------------------------------------------------------------
        private static double MinJanuaryTempModifier(AnnualClimate_Monthly weather, ISpecies species)
        // Is the January mean temperature greater than the species specified minimum?
        {
        
            int speciesMinimum = SpeciesData.MinJanTemp[species];
            
            if (weather.MonthlyTemp[0] < speciesMinimum)
                return 0.0;
            else
                return 1.0;
        }
    }
}
