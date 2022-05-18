using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using Rhino.Geometry.Intersect;
using System.Linq;

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
        public static Point3d m_SrfPt;
        public static List<Point3d> m_SrfPts;
        public static double m_minU;
        public static double m_maxU;
        public static double m_minV;
        public static double m_maxV;

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            //pManager.AddNumberParameter("d", "d", "d", GH_ParamAccess.list);
            //pManager.AddNumberParameter("side", "side", "side", GH_ParamAccess.item);
            pManager.AddNumberParameter("u", "u", "u", GH_ParamAccess.item, 0);
            pManager.AddNumberParameter("v", "v", "v", GH_ParamAccess.item, 0);
            pManager.AddNumberParameter("u_scale", "u_scale", "u_scale", GH_ParamAccess.item, 0.5);
            pManager.AddNumberParameter("v_scale", "v_scale", "v_scale", GH_ParamAccess.item, 0.5);
            pManager[2].Optional = true;
            pManager[3].Optional = true;
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
            pManager.AddGenericParameter("kato", "kato", "kato", GH_ParamAccess.item);
            pManager.AddGenericParameter("kato_pointlist", "kato_pointlist", "kato_pointlist", GH_ParamAccess.list);

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
            int Side = 0;

            //DA.GetDataList(0, Dvalues);
            //DA.GetData(1,ref Side);

            double u = 0, v = 0, u_scl = 0.1, v_scl = 0.1;
            DA.GetData(0, ref u);
            DA.GetData(1, ref v);
            DA.GetData(2, ref u_scl);
            DA.GetData(3, ref v_scl);
            if (u_scl == 0 || v_scl == 0) return;

            double Kappa = 0;

            double landa = ComputeSideBLFunction(Dvalues, Side);
            if (Side != 0)
                Kappa = ComputeCornerBLFunction(Dvalues, Side - 1, Side);
            else
                Kappa = ComputeCornerBLFunction(Dvalues, Dvalues.Count - 1, Side);

            double Nu = ComputeSPsideBLFunction(Dvalues, Side);

            ComputeDomainPolygon();

            m_SrfPt = Kato_Suv(u, v);
            ComputeSrfPts(u_scl, v_scl);
            var num_u = 20;
            var num_v = 20;
            ComputeSurfacePoints(num_u, num_v);

            DA.SetDataList(0, m_Curves);
            DA.SetDataList(1, m_DPolygonLines);
            DA.SetData(2, landa);
            DA.SetData(3, Kappa);
            DA.SetData(4, Nu);
            DA.SetData(5, m_SrfPt);
            DA.SetDataList(6, m_SrfPts);

        }


        private void ComputeSurfacePoints(double u_div, double v_div)
        {
            ComputeMinMaxUVs();
            var polygon = new List<Point3d>(m_DomainPolygon);
            Point3d first = new Point3d(polygon[0]);
            polygon.Add(first);
            var poly = Curve.CreateInterpolatedCurve(polygon, 1);

            var inc_U = (double)(m_maxU - m_minU) / u_div;
            var inc_V = (double)(m_maxV - m_minV) / v_div;
            var intlist = new List<int>();
            m_SrfPts = new List<Point3d>();
            for (double i = m_minU; i <= m_maxU; i += inc_U)
            {
                var counterx = 0;
                var Ln = new Line(new Point3d(i, m_minV, 0), new Point3d(i, m_maxV, 0));
                var sr = Intersection.CurveLine(poly, Ln, 0.0001, 0.001);
                for (double j = m_minV; j <= m_maxV; j += inc_V)
                {
                    var state = poly.Contains(new Point3d(i, j, 0), Plane.WorldXY, 0.0000000001);
                    counterx++;
                    if (state != PointContainment.Outside)
                    {
                        //if (m_DomainPolygon.Contains(new Point3d(i, j, 0)))
                        //    //m_SrfPts.Add(m_Curves[m_DomainPolygon.IndexOf(new Point3d(i, j, 0))].PointAtStart);
                        //else
                        m_SrfPts.Add(Kato_Suv(i, j));
                    }
                    else
                    {
                        var new_pt = new Point3d(i, j, 0);
                        var dist_list = new List<double>();
                        var edge_point = new List<Point3d>();

                        List<Point3d> pts = new List<Point3d>();

                        for (int k = 0; k < sr.Count; k++)
                        {
                            edge_point.Add(sr[k].PointA);
                            dist_list.Add(new_pt.DistanceTo(sr[k].PointA));
                        }

                        var min_val = dist_list.Min();
                        var idx = dist_list.IndexOf(min_val);


                        //if (m_DomainPolygon.Contains(edge_point[idx]))
                        //{
                        //    var value = m_Curves[m_DomainPolygon.IndexOf(edge_point[idx])].PointAtStart;
                        //    //var value =  Kato_Suv(edge_point[idx].X, edge_point[idx].Y);
                        //    m_SrfPts.Add(Kato_Suv(edge_point[idx].X, edge_point[idx].Y));

                        //    var sd = 1;
                        //}
                        ////else
                        m_SrfPts.Add(Kato_Suv(edge_point[idx].X, edge_point[idx].Y));
                    }

                }
                intlist.Add(counterx);
            }
        }








        private void ComputeSrfPts(double u_scl, double v_scl)
        {
            int s;
            double N;
            Point3d Pt;

            ComputeMinMaxUVs();

            m_SrfPts = new List<Point3d>();
            for (double i = m_minU; i <= m_maxU; i += u_scl)
            {
                for (double j = m_minV; j <= m_maxV; j += v_scl)
                {
                    if (!IsInDomainPolygon(i, j))
                        continue;

                    if (Kato_Suv(i, j) != null)
                        m_SrfPts.Add(Kato_Suv(i, j));
                    //j += 0.05;
                }
            }

            N = 10.0;
            for (int i = 0; i < m_DomainPolygon.Count; i++)
            {
                s = (i + 1) % m_DomainPolygon.Count;

                for (double k = 0; k < N; k = k + 1.0)
                {
                    Pt = m_DomainPolygon[i] * (N - k) / N + m_DomainPolygon[s] * k / N;
                    var val = Kato_Suv(Pt.X, Pt.Y);

                    if (val.IsValid)
                        m_SrfPts.Add(Kato_Suv(Pt.X, Pt.Y));
                }
            }
        }

        private void ComputeMinMaxUVs()
        {
            m_minU = 9E9;
            m_maxU = -9E9;
            m_minV = 9E9;
            m_maxV = -9E9;

            for (int i = 0; i < m_DomainPolygon.Count; i++)
            {
                if (m_DomainPolygon[i].X < m_minU)
                    m_minU = m_DomainPolygon[i].X;
                if (m_DomainPolygon[i].Y < m_minV)
                    m_minV = m_DomainPolygon[i].Y;

                if (m_DomainPolygon[i].X > m_maxU)
                    m_maxU = m_DomainPolygon[i].X;
                if (m_DomainPolygon[i].Y > m_maxV)
                    m_maxV = m_DomainPolygon[i].Y;
            }
        }



        private bool IsInDomainPolygon(double u, double v)
        {
            List<Vector3d> Vecs = new List<Vector3d>();
            Vector3d Vec;
            int j;
            double TotAng;




            for (int i = 0; i < m_DomainPolygon.Count; i++)
            {
                Vec = new Vector3d(u - m_DomainPolygon[i].X, v - m_DomainPolygon[i].Y, 0);
                Vecs.Add(Vec);
            }

            TotAng = 0;
            for (int i = 0; i < Vecs.Count; i++)
            {
                j = (i + 1) % Vecs.Count;
                TotAng = TotAng + Vector3d.VectorAngle(Vecs[i], Vecs[j]);
            }

            if (Math.Abs(TotAng - 2 * Math.PI) < 9E-9)
                return true;

            return false;
        }

        private Point3d Kato_Suv(double u, double v)
        {
            List<(double, double)> si_di = ComputeRadialDistanceFunctionforKato(u, v, m_Curves.Count);

            List<double> d_i = new List<double>();
            for (int i = 0; i < si_di.Count; i++)
            {
                d_i.Add(si_di[i].Item2);
            }

            Point3d r_sum = new Point3d();
            List<double> Value = ComputeSPsideBLFunction1(d_i);
            for (int i = 0; i < m_Curves.Count; i++)
            {
                double domainlengt = (m_DomainPolygon[IndexWrapper(i, m_Curves.Count)] - m_DomainPolygon[IndexWrapper(i + 1, m_Curves.Count)]).Length;
                double curveparameter = (m_Curves[i].Domain.Length * si_di[i].Item1) / domainlengt;

                //Vector3d crossproduct = Vector3d.CrossProduct(m_Curves[i].TangentAt(curveparameter), m_Curves[i].CurvatureAt(curveparameter));
                //Point3d r = m_Curves[i].PointAt(curveparameter) + (si_di[i].Item2 * crossproduct);

                Point3d r = m_Curves[i].PointAt(curveparameter) + (si_di[i].Item2 * m_Curves[i].TangentAt(curveparameter));
                //double blendingfunctionvalue = ComputeSPsideBLFunction(d_i, i);
                r_sum += r * Value[i];
            }
            return r_sum;
        }

        private List<(double, double)> ComputeRadialDistanceFunctionforKato(double u, double v, int n)
        {
            List<(double, double)> output = new List<(double, double)>();
            int i_before = 0, i_after = 0, i_afterafter = 0;

            for (int i = 0; i < n; i++)
            {

                i_before = IndexWrapper((i - 1), n);
                i_after = IndexWrapper((i + 1), n);
                i_afterafter = IndexWrapper((i + 2), n);

                Line i_beforeline = new Line(m_DomainPolygon[i], m_DomainPolygon[i_before]);
                Line i_afterline = new Line(m_DomainPolygon[i_after], m_DomainPolygon[i_afterafter]);

                bool A = Intersection.LineLine(i_beforeline, i_afterline, out double a, out double b, 0.001, false);
                Point3d c_i = i_beforeline.PointAt(a);
                Point3d c_i_0 = i_afterline.PointAt(b);


                Line null_line = new Line(new Point3d(u, v, 0), c_i);
                Line null_line2 = new Line(m_DomainPolygon[i], m_DomainPolygon[i_after]);
                bool B = Intersection.LineLine(null_line, null_line2, out double a1, out double b2, 0.001, true);
                Point3d e_i = null_line.PointAt(a1);
                Point3d e_i0 = null_line2.PointAt(b2);

                if (!A || !B || (c_i - c_i_0).Length > 0.000001 || (e_i - e_i0).Length > 0.000001) AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Radial distance function domain line Intersection problem!!");

                //output.Add(((new Point3d(u, v, 0) - e_i).Length, (e_i - m_DomainPolygon[i]).Length));
                output.Add(((e_i - m_DomainPolygon[i]).Length, (new Point3d(u, v, 0) - e_i).Length));
            }
            return output;
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
            for (int i = 0; i < m_DomainPolygon.Count - 1; i++)
                m_DPolygonLines.Add(new Line(m_DomainPolygon[i], m_DomainPolygon[i + 1]));
            m_DPolygonLines.Add(new Line(m_DomainPolygon[m_DomainPolygon.Count - 1], m_DomainPolygon[0]));
        }

        private static double ComputeSideBLFunction(List<double> d, int Side)
        {
            double Landa = 0;

            double Numerator = 0;
            double denominator = 0;
            List<int> Numerator1 = new List<int>();
            List<int> Numerator2 = new List<int>();



            Numerator1.Add(Side);
            Numerator2.Add(Side);



            if (Side == 0)
                Numerator1.Add(d.Count - 1);
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

        private static double ComputeCornerBLFunction(List<double> d, int Side1, int Side2) //Insert Side in order from 0 to n
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

        private static List<double> ComputeSPsideBLFunction1(List<double> d) //Insert Side in order from 0 to n
        {
            List<double> Mus = new List<double>();
            for (int i = 0; i < d.Count; i++)
            {
                double Nu = 0;

                double Numerator = 0;
                double denominator = 0;
                List<int> Numerator1 = new List<int>();

                Numerator1.Add(i);

                Numerator = ProductFunction(d, Numerator1);

                for (int j = 0; j < d.Count; j++)
                {
                    List<int> Denominator1 = new List<int>();
                    Denominator1.Add(j);

                    double ProVal = ProductFunction(d, Denominator1);
                    denominator += ProVal;

                }

                Nu = Numerator / denominator;

                if (double.IsNaN(Nu))
                    Nu = 0;
                Mus.Add(Nu);

            }

            for (int i = 0; i < Mus.Count; i++)
            {
                var side = i;
                var Preside = i - 1;

                if (i == 0)
                    Preside = d.Count - 1;

                if (d[side] < 0.0001 && d[Preside] < 0.0001)
                {
                    Mus[side] = 1;
                    Mus[Preside] = 0;

                }

            }


            return Mus;
        }
        private static double ProductFunction(List<double> d, List<int> n)
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
                SqreList.Add(Math.Pow(d[i], 2));
            }

            double Product = 1.0;
            for (int i = 0; i < SqreList.Count; i++)
            {
                Product *= SqreList[i];
            }

            Result = Product;

            return Result;

        }

        /// <summary>
        /// List Modulo Operator for all integer index value
        /// </summary>
        /// <param name="index"></param>
        /// <param name="list_count"></param>
        /// <returns></returns>
        private int IndexWrapper(int index, int list_count)
        {
            return ((index % list_count) + list_count) % list_count;
        }
        private static double ToRad(double Deg)
        {
            return Math.PI * Deg / 180.0;
        }

        private static double ToDeg(double Rad)
        {
            return Rad * 180.0 / Math.PI;
        }
        public override void AddRuntimeMessage(GH_RuntimeMessageLevel level, string text)
        {
            base.AddRuntimeMessage(level, text);
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