//  Copyright 2007-2010 Portland State University, University of Wisconsin-Madison
//  Author: Robert Scheller, Melissa Lucash

using Edu.Wisc.Forest.Flel.Util;
using System;
using System.Threading;
using Landis.Core;
using Landis.SpatialModeling;


namespace Landis.Extension.Succession.NECN_Hydro
{

    public enum LayerName {Leaf, FineRoot, Wood, CoarseRoot, Metabolic, Structural, SOM1, SOM2, SOM3, Other};
    public enum LayerType {Surface, Soil, Other} 

    /// <summary>
    /// A Century soil model carbon and nitrogen pool.
    /// </summary>
    public class Layer
    {
        private LayerName name;
        private LayerType type;
        private double carbon;
        private double nitrogen;
        private double decayValue;
        private double fractionLignin;
        private double netMineralization;
        private double grossMineralization;


        //---------------------------------------------------------------------
        public Layer(LayerName name, LayerType type)
        {
            this.name = name;
            this.type = type;
            this.carbon = 0.0;
            this.nitrogen = 0.0;

            this.decayValue = 0.0;
            this.fractionLignin = 0.0;

            this.netMineralization = 0.0;
            this.grossMineralization = 0.0;

        }
        //---------------------------------------------------------------------
        /// <summary>
        /// Name
        /// </summary>
        public LayerName Name { get; set; }
        
        /// <summary>
        /// Provides an index to LitterTypeTable
        /// </summary>
        public LayerType Type { get; set; }
        
        /// <summary>
        /// Carbon
        /// </summary>
        public double Carbon { get; set; }
        
        /// <summary>
        /// Nitrogen
        /// </summary>
        public double Nitrogen { get; set; }
        
        /// <summary>
        /// Pool decay rate.
        /// </summary>
        public double DecayValue { get; set; }
        
        //---------------------------------------------------------------------
        /// <summary>
        /// Pool Carbon:Nitrogen Ratio
        /// </summary>
        public double FractionLignin { get; set; }
        
        /// <summary>
        /// Net Mineralization
        /// </summary>
        public double NetMineralization { get; set; }
        
        /// <summary>
        /// Gross Mineralization
        /// </summary>
        public double GrossMineralization { get; set; }
        
        public Layer Clone()
        {
            Layer newLayer = new Layer(this.Name, this.Type);

            newLayer.carbon = this.carbon;
            newLayer.nitrogen = this.nitrogen ;

            newLayer.decayValue = this.decayValue ;
            newLayer.fractionLignin = this.fractionLignin ;

            newLayer.netMineralization = this.netMineralization ;
            newLayer.grossMineralization = this.grossMineralization ;

            return newLayer;
        }

        // --------------------------------------------------
        public void DecomposeStructural(ActiveSite site)
        {
            if (this.Carbon > 0.0000001)
            {

                double anerb = SiteVars.AnaerobicEffect[site];

                if (this.Type == LayerType.Surface) anerb = 1.0; // No anaerobic effect on surface material

                //Compute total C flow out of structural in layer
                double totalCFlow = System.Math.Min(this.Carbon, OtherData.MaxStructuralC)
                                * SiteVars.DecayFactor[site]
                                * OtherData.LitterParameters[(int)this.Type].DecayRateStrucC
                                * anerb
                                * System.Math.Exp(-1.0 * OtherData.LigninDecayEffect * this.FractionLignin)
                                * OtherData.MonthAdjust;

                //Decompose structural into som1 and som2 with CO2 loss.
                this.DecomposeLignin(totalCFlow, site);
            }
        }
        // --------------------------------------------------
        public void DecomposeLignin(double totalCFlow, ActiveSite site)
        // Originally from declig.f for decomposition of compartment lignin
        {
            //double carbonToSOM1;    //Net C flow to SOM1
            //double carbonToSOM2;    //Net C flow to SOM2
            double litterC = this.Carbon; 
            double ratioCN = litterC / this.Nitrogen;


            //See if Layer can decompose to SOM1.
            //If it can decompose to SOM1, it will also go to SOM2.
            //If it can't decompose to SOM1, it can't decompose at all.

            //If Wood Object can decompose
            if (this.DecomposePossible(ratioCN, SiteVars.MineralN[site]))
            {
                // Decompose Wood Object to SOM2
                // -----------------------
                // Gross C flow to som2
                double carbonToSOM2 = totalCFlow * this.FractionLignin;

                //MicrobialRespiration associated with decomposition to som2
                double co2loss = carbonToSOM2 * OtherData.LigninRespirationRate;

                this.Respiration(co2loss, site);

                //Net C flow to SOM2
                double netCFlow = carbonToSOM2 - co2loss;

                // Partition and schedule C flows 
                this.TransferCarbon(SiteVars.SOM2[site], netCFlow);
                this.TransferNitrogen(SiteVars.SOM2[site], netCFlow, litterC, ratioCN, site);

                // ----------------------------------------------
                // Decompose Wood Object to SOM1
                // Gross C flow to som1

                double carbonToSOM1 = totalCFlow - carbonToSOM2 - co2loss;

                //MicrobialRespiration associated with decomposition to som1
                if(this.Type == LayerType.Surface)
                    co2loss = carbonToSOM1 * OtherData.StructuralToCO2Surface; 
                else
                    co2loss = carbonToSOM1 * OtherData.StructuralToCO2Soil; 

                this.Respiration(co2loss, site);

                //Net C flow to SOM1
                carbonToSOM1 -= co2loss;

                if(this.Type == LayerType.Surface)
                {
                    this.TransferCarbon(SiteVars.SOM1surface[site], carbonToSOM1);
                    this.TransferNitrogen(SiteVars.SOM1surface[site], carbonToSOM1, litterC, ratioCN, site);
                }
                else
                {
                    this.TransferCarbon(SiteVars.SOM1soil[site], carbonToSOM1);
                    this.TransferNitrogen(SiteVars.SOM1soil[site], carbonToSOM1, litterC, ratioCN, site);
                }
            }
        }

        //---------------------------------------------------------------------
        public void DecomposeMetabolic(ActiveSite site)
        {
            double litterC = this.Carbon;
            double anerb = SiteVars.AnaerobicEffect[site];

            if (litterC > 0.0000001)
            {
              // Determine C/N ratios for flows to SOM1
                double ratioCNtoSOM1 = 0.0;
                double co2loss = 0.0;

                // Compute ratios for surface  metabolic residue
                if (this.Type == LayerType.Surface)
                    ratioCNtoSOM1 = Layer.AbovegroundDecompositionRatio(this.Nitrogen, litterC);

                //Compute ratios for soil metabolic residue
                else
                    ratioCNtoSOM1 = Layer.BelowgroundDecompositionRatio(site,
                                        OtherData.MinCNenterSOM1,
                                        OtherData.MaxCNenterSOM1,
                                        OtherData.MinContentN_SOM1);

                //Compute total C flow out of metabolic layer
                double totalCFlow = litterC
                                * SiteVars.DecayFactor[site]
                                * OtherData.LitterParameters[(int) this.Type].DecayRateMetabolicC
                                * OtherData.MonthAdjust;

                //Added impact of soil anerobic conditions
                if (this.Type == LayerType.Soil) totalCFlow *= anerb;

                //Make sure metabolic C does not go negative.
                if (totalCFlow > litterC)
                    totalCFlow = litterC;

                //If decomposition can occur,
                if (this.DecomposePossible(ratioCNtoSOM1, SiteVars.MineralN[site]))
                {
                    //CO2 loss
                    if (this.Type == LayerType.Surface)
                        co2loss = totalCFlow * OtherData.MetabolicToCO2Surface;
                    else
                        co2loss = totalCFlow * OtherData.MetabolicToCO2Soil;

                    this.Respiration(co2loss, site);

                    //Decompose metabolic into som1
                    double netCFlow = totalCFlow - co2loss;

                    if (netCFlow > litterC)
                        PlugIn.ModelCore.UI.WriteLine("   ERROR:  Decompose Metabolic:  netCFlow={0:0.000} > layer.Carbon={0:0.000}.", netCFlow, this.Carbon);

                    // -- CARBON AND NITROGEN ---------------------------
                    // Partition and schedule C flows
                    // Compute and schedule N flows and update mineralization accumulators.
                    if((int) this.Type == (int) LayerType.Surface)
                    {
                        this.TransferCarbon(SiteVars.SOM1surface[site], netCFlow);
                        this.TransferNitrogen(SiteVars.SOM1surface[site], netCFlow, litterC, ratioCNtoSOM1, site);
                    }
                    else
                    {
                        this.TransferCarbon(SiteVars.SOM1soil[site], netCFlow);
                        this.TransferNitrogen(SiteVars.SOM1soil[site], netCFlow, litterC, ratioCNtoSOM1, site);
                    }

                }
            }
        }

        //---------------------------------------------------------------------
        public void TransferCarbon(Layer destination, double netCFlow)
        {

            if (netCFlow > this.Carbon)
            {
                netCFlow = this.Carbon;
            }
                
            this.Carbon -= netCFlow;
            destination.Carbon += netCFlow;
        }

        public void TransferNitrogen(Layer destination, double CFlow, double totalC, double ratioCNtoDestination, ActiveSite site)
        {
            double mineralNFlow = 0.0;

            //...N flow is proportional to C flow.
            double NFlow = this.Nitrogen * CFlow / totalC;

            //...This was added to avoid a 0/0 error on the pc.
            if (CFlow <= 0.0 || NFlow <= 0.0)
            {
                return;
            }

            //...If C/N of Box A > C/N of new material entering Box B
            if ((CFlow / NFlow) > ratioCNtoDestination)
            {
                //...IMMOBILIZATION occurs.
                //...Compute the amount of N immobilized.
                //     since  ratioCNtoDestination = netCFlow / (Nflow + immobileN),
                //     where immobileN is the extra N needed from the mineral pool
                double immobileN = (CFlow / ratioCNtoDestination) - NFlow;
                
                this.Nitrogen -= NFlow;
                destination.Nitrogen += NFlow;

                //Don't allow mineral N to go to zero or negative.- ML
                if (immobileN > SiteVars.MineralN[site])
                {
                    immobileN = SiteVars.MineralN[site] - 0.01; //leave some small amount of mineral N
                }

                SiteVars.MineralN[site] -= immobileN;
                destination.Nitrogen += immobileN;

                //...Return mineralization value.
                mineralNFlow = -1 * immobileN;
            }
            else
            {
                //...MINERALIZATION occurs
                //...Schedule flow from Box A to Box B
                double mineralizedN = (CFlow / ratioCNtoDestination);

                this.Nitrogen -= mineralizedN;
                destination.Nitrogen += mineralizedN;

                //...Schedule flow from Box A to mineral pool

                mineralNFlow = NFlow - mineralizedN;

                this.Nitrogen -= mineralNFlow;

                SiteVars.MineralN[site] += mineralNFlow;
            }

            if (mineralNFlow > 0)
            {
                SiteVars.GrossMineralization[site] += mineralNFlow;
            }

            //...Net mineralization
            this.NetMineralization += mineralNFlow;
        }

        public void Respiration(double co2loss, ActiveSite site)
        {
            // Copyright 1993 Colorado State University
            // All Rights Reserved
            // Compute flows associated with microbial respiration.

            // Input:
            //  co2loss = CO2 loss associated with decomposition
            //  Box A.  For components with only 1 layer, tcstva will be dimensioned (1).

            //  Transput:
            //c         carbonSourceSink = C source/sink
            //c         grossMineralization = gross mineralization
            //c         netMineralization = net mineralization for layer N

            //c...Mineralization associated with respiration is proportional to the N fraction.
            double mineralNFlow = co2loss * this.Nitrogen / this.Carbon; 

            if(mineralNFlow > this.Nitrogen)
            {
                mineralNFlow = this.Nitrogen;
                co2loss = this.Carbon;
            }

            this.TransferCarbon(SiteVars.SourceSink[site], co2loss);

            //Add lost CO2 to monthly heterotrophic respiration
            SiteVars.MonthlyResp[site][Century.Month] += co2loss;

            this.Nitrogen -= mineralNFlow;
            SiteVars.MineralN[site] += mineralNFlow;

            //c...Update gross mineralization
            if (mineralNFlow > 0)
                SiteVars.GrossMineralization[site] += mineralNFlow;

            //c...Update net mineralization
            this.NetMineralization += mineralNFlow;
        }

        public bool DecomposePossible(double ratioCNnew, double mineralN)
        {

            //Copyright 1993 Colorado State University
            //All Rights Reserved
            //logical function candec(nelem,aminrl,tca,elstva,nlr,lyr,rcenew)
            //c...Determine if decomposition can occur.

            bool canDecompose = true;

            //c...If there is no available mineral N
            
            // Compare the C/N of new material to the C/N of the layer if C/N of
            // the layer > C/N of new material
            if (mineralN < 0.0000001 && (this.Carbon / this.Nitrogen > ratioCNnew))
            {
                // Immobilization is necessary and the stuff in Box A can't
                // decompose to Box B.
                canDecompose = false;
            }

            // If there is some available mineral N, decomposition can
            // proceed even if mineral N is driven negative in
            // the next time step.

            return canDecompose;
        }

        public void AdjustLignin(double inputC, double inputFracLignin)
        {
            // Copyright 1993 Colorado State University
            // All Rights Reserved

            //c...Adjust the fraction of lignin in structural C when new material
            //c...  is added.

            //c    oldc  = grams C in structural before new material is added
            //c    frnew = fraction of lignin in new structural material
            //c    addc  = grams structural C being added

            //c    fractl comes in as fraction of lignin in structural before new
            //c           material is added; goes out as fraction of lignin in
            //c           structural with old and new combined.

            //c...oldlig  = grams of lignin in existing residue
            double oldlig = this.FractionLignin * this.Carbon;//totalC;

            //c...newlig = grams of lignin in new residue
            double newlig = inputFracLignin * inputC;

            //c...Compute lignin fraction in combined residue
            double newFraction = (oldlig + newlig) / (this.Carbon + inputC);

            this.FractionLignin = newFraction;
        }

        public void AdjustDecayRate(double inputC, double inputDecayRate)
        {
            //c...oldlig  = grams of lignin in existing residue
            double oldDecayRate = this.DecayValue * this.Carbon;

            //c...newlig = grams of lignin in new residue
            double newDecayRate = inputDecayRate * inputC;

            //c...Compute decay rate in combined residue
            this.DecayValue = (oldDecayRate + newDecayRate) / (inputC + this.Carbon);
        }


        public static double BelowgroundDecompositionRatio(ActiveSite site, double minCNenter, double maxCNenter, double minContentN)
        {
            //Originally from bgdrat.f
            //BelowGround Decomposition RATio computation.
            double bgdrat = 0.0;

            //Determine ratio of C/N of new material entering 'Box B'.
            //Ratio depends on available N

            double mineralN = SiteVars.MineralN[site];

            if (mineralN <= 0.0)
                bgdrat = maxCNenter;  // Set ratio to maximum allowed (HIGHEST carbon, LOWEST nitrogen)
            else if (mineralN > minContentN)
                bgdrat = minCNenter;  //Set ratio to minimum allowed
            else
                bgdrat = (1.0 - (mineralN / minContentN)) * (maxCNenter - minCNenter)
                    + minCNenter;

            return bgdrat;
        }

        public static double AbovegroundDecompositionRatio(double abovegroundN, double abovegroundC)
        {       //Originally from agdrat.f.

            double Ncontent, agdrat;
            double biomassConversion = 2.0;
            // cemicb = slope of the regression line for C/N of som1

            double cemicb = (OtherData.MinCNSurfMicrobes - OtherData.MaxCNSurfMicrobes) / OtherData.MinNContentCNSurfMicrobes;


            //The C/E ratios for structural and wood can be computed once;
            //they then remain fixed throughout the run.  The ratios for
            //metabolic and som1 may vary and must be recomputed each time step

            if ((abovegroundC * biomassConversion) <= 0.00000000001)  Ncontent = 0.0;
            else  Ncontent = abovegroundN / (abovegroundC * biomassConversion);

            //tca is multiplied by biomassConversion to give biomass

            if (Ncontent > OtherData.MinNContentCNSurfMicrobes)
                agdrat = OtherData.MinCNSurfMicrobes;
            else
                agdrat = OtherData.MaxCNSurfMicrobes + Ncontent * cemicb;

            return agdrat;
        }

        //---------------------------------------------------------------------
        /// <summary>
        /// Reduces the pool's biomass by a specified percentage.
        /// </summary>
        public void ReduceMass(double percentageLost)
        {
            if (percentageLost < 0.0 || percentageLost > 1.0)
                throw new ArgumentException("Percentage must be between 0% and 100%");

            this.Carbon   = this.Carbon * (1.0 - percentageLost);
            this.Nitrogen = this.Nitrogen * (1.0 - percentageLost);
        }

    }
}
