using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;

namespace Grasshopper2
{
    public class Transfinite : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the Transfinite class.
        /// </summary>
        public Transfinite()
          : base("Transfinite", "Transfinite",
              "Transfinite",
              "Transfinite", "Transfinite")
        {
        }

        private static List<Curve> m_Curves;
        private static List<Point3d> m_DomainPolygon;
        private static List<Line> m_DPolygonLines;
        public static List<double> m_DisToEdges;

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddNumberParameter("d", "d", "d", GH_ParamAccess.list);
            pManager.AddNumberParameter("side", "side", "side", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddCurveParameter("Curves", "Curves", "Curves", GH_ParamAccess.list);
            pManager.AddLineParameter("DomainPolygon", "DomainPolygon", "DomainPolygon", GH_ParamAccess.list);
            pManager.AddNumberParameter("landa", "landa", "landa", GH_ParamAccess.item);
            pManager.AddNumberParameter("kappa", "kappa", "kappa", GH_ParamAccess.item);
            pManager.AddNumberParameter("nu", "nu", "nu", GH_ParamAccess.item);

        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {

            List<Point3d> Pts1 = new List<Point3d>() { new Point3d(0, 0, 0), new Point3d(3, 0, 1), new Point3d(6, 0, -1), new Point3d(10, -1, 0) };
            List<Point3d> Pts2 = new List<Point3d>() { new Point3d(10, -1, 0), new Point3d(12, 3, -1), new Point3d(14, 7, 1), new Point3d(15, 10, 0) };
            List<Point3d> Pts3 = new List<Point3d>() { new Point3d(15, 10, 0), new Point3d(13, 13, -1), new Point3d(12, 16, -1), new Point3d(10, 20, 0) };
            List<Point3d> Pts4 = new List<Point3d>() { new Point3d(10, 20, 0), new Point3d(7, 17, 1), new Point3d(3, 14, -1), new Point3d(0.5, 11, 0) };
            List<Point3d> Pts5 = new List<Point3d>() { new Point3d(0.5, 11, 0), new Point3d(0, 8, 1), new Point3d(0, 4, -1), new Point3d(0, 0, 0) };

            Curve Crv1 = Curve.CreateInterpolatedCurve(Pts1, 3);
            Curve Crv2 = Curve.CreateInterpolatedCurve(Pts2, 3);
            Curve Crv3 = Curve.CreateInterpolatedCurve(Pts3, 3);
            Curve Crv4 = Curve.CreateInterpolatedCurve(Pts4, 3);
            Curve Crv5 = Curve.CreateInterpolatedCurve(Pts5, 3);
            m_Curves = new List<Curve>();
            m_Curves.Add(Crv1);
            m_Curves.Add(Crv2);
            m_Curves.Add(Crv3);
            m_Curves.Add(Crv4);
            m_Curves.Add(Crv5);

            List<double> Dvalues = new List<double>();
            int Side=0;
            
            DA.GetDataList(0, Dvalues);
            DA.GetData(1,ref Side);

            double Kappa = 0;

            double landa = ComputeSideBLFunction(Dvalues, Side);
            if (Side != 0)
            {
                Kappa = ComputeCornerBLFunction(Dvalues, Side - 1, Side);
            }

            else
            {
                Kappa = ComputeCornerBLFunction(Dvalues, Dvalues.Count - 1, Side);
            } 
                
           
            double Nu = ComputeSPsideBLFunction(Dvalues, Side);

            ComputeDomainPolygon();

            DA.SetDataList(0, m_Curves);
            DA.SetDataList(1, m_DPolygonLines);

            DA.SetData(2, landa);
            DA.SetData(3, Kappa);
            DA.SetData(4, Nu);
        }

        private static void ComputeDomainPolygon()
        {
            double TotL = 0, L;
            List<double> Ls = new List<double>();
            for (int i = 0; i < m_Curves.Count; i++)
            {
                L = m_Curves[i].GetLength();
                TotL = TotL + L;
                Ls.Add(L);
            }            
                

            List<double> Alphas = new List<double>();
            Alphas.Add(0.0);
            for (int i = 1; i < m_Curves.Count; i++)
            {
                L = 0;
                for (int j = 0; j <= i - 1; j++)
                    L = L + Ls[j];
                Alphas.Add(2 * Math.PI * L / TotL);
            }
            
            m_DomainPolygon = new List<Point3d>();
            Point3d Pt;
            for (int i = 0; i < Alphas.Count; i++)
            {
                Pt = new Point3d(Math.Cos(Alphas[i]), Math.Sin(Alphas[i]), 0);
                m_DomainPolygon.Add(Pt);
            }

            m_DPolygonLines = new List<Line>();
            for (int i = 0; i < m_DomainPolygon.Count-1; i++)
                m_DPolygonLines.Add(new Line(m_DomainPolygon[i], m_DomainPolygon[i+1]));
            m_DPolygonLines.Add(new Line(m_DomainPolygon[m_DomainPolygon.Count-1], m_DomainPolygon[0]));
        }

        private static double ComputeSideBLFunction(List<double> d,int Side)
        {
            double Landa = 0;

            double Numerator = 0;
            double denominator = 0;
            List<int> Numerator1 = new List<int>();
            List<int> Numerator2 = new List<int>();

            

            Numerator1.Add(Side);
            Numerator2.Add(Side);



            if (Side == 0)
                Numerator1.Add(d.Count-1);
            else
                Numerator1.Add(Side - 1);

            if (Side == d.Count - 1)
                Numerator2.Add(0);
            else
                Numerator2.Add(Side + 1);

            Numerator = ProductFunction(d, Numerator1) + ProductFunction(d, Numerator2);

            for (int i = 0; i < d.Count; i++)
            {
                List<int> Denominator1 = new List<int>();
                Denominator1.Add(i);

                if (i == 0)
                {
                    Denominator1.Add(d.Count - 1);
                }    
                else
                {
                    Denominator1.Add(i - 1);
                }
                    

                double ProVal = ProductFunction(d, Denominator1);
                denominator += ProVal;

            }

            Landa = Numerator / denominator;

            return Landa; 
        }

        private static double ComputeCornerBLFunction(List<double> d, int Side1 , int Side2) //Insert Side in order from 0 to n
        {
            double Kapa = 0;

            double Numerator = 0;
            double denominator = 0;
            List<int> Numerator1 = new List<int>();

            Numerator1.Add(Side1);
            Numerator1.Add(Side2);


            Numerator = ProductFunction(d, Numerator1);

            for (int i = 0; i < d.Count; i++)
            {
                List<int> Denominator1 = new List<int>();
                Denominator1.Add(i);

                if (i == 0)
                    Denominator1.Add(d.Count - 1);
                else
                    Denominator1.Add(i - 1);

                double ProVal = ProductFunction(d, Denominator1);
                denominator += ProVal;

            }

            Kapa = Numerator / denominator;

            return Kapa;
        }

        private static double ComputeSPsideBLFunction(List<double> d, int Side) //Insert Side in order from 0 to n
        {
            double Nu = 0;

            double Numerator = 0;
            double denominator = 0;
            List<int> Numerator1 = new List<int>();

            Numerator1.Add(Side);

            Numerator = ProductFunction(d, Numerator1);

            for (int i = 0; i < d.Count; i++)
            {
                List<int> Denominator1 = new List<int>();
                Denominator1.Add(i);

                double ProVal = ProductFunction(d, Denominator1);
                denominator += ProVal;

            }

            Nu = Numerator / denominator;

            return Nu;
        }

        private static double ProductFunction(List<double> d,List<int> n)
        {
            double Result = 0;
            List<double> SqreList = new List<double>();

            for (int i = 0; i < d.Count; i++)
            {
                bool eliminated = false;

                for (int j = 0; j < n.Count; j++)
                {
                    if (i == n[j])
                    {
                        eliminated = true;
                        break;
                    }   
                }

                if (eliminated)
                    continue;
                SqreList.Add(Math.Pow(d[i],2));
            }

            double Product = 1.0;
            for (int i = 0; i < SqreList.Count; i++)
            {
                Product *= SqreList[i]; 
            }

            Result = Product;

            return Result;

        }

        private static double ToRad(double Deg)
        {
            return Math.PI * Deg / 180.0;
        }

        private static double ToDeg(double Rad)
        {
            return Rad * 180.0 / Math.PI;
        }

        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                // return Resources.IconForThisComponent;
                return null;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("35212335-292f-402d-8f09-1e0290d10f4a"); }
        }
    }
}