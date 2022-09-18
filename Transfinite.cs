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
          : base("Transfinite", "Nickname",
            "Description",
            "Transfinite", "Subcategory")
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
        private static List<Brep> m_Patches;
        private static List<double> m_PlanarityErrors;
        private static double m_PlnAverageError;
        public static List<Point3d> m_uv;

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            //pManager.AddNumberParameter("d", "d", "d", GH_ParamAccess.list);
            //pManager.AddNumberParameter("side", "side", "side", GH_ParamAccess.item);
            pManager.AddNumberParameter("u", "u", "u", GH_ParamAccess.item, 0);
            pManager.AddNumberParameter("v", "v", "v", GH_ParamAccess.item, 0);
            pManager.AddIntegerParameter("u_scale", "u_scale", "u_scale", GH_ParamAccess.item, 40);
            pManager.AddIntegerParameter("v_scale", "v_scale", "v_scale", GH_ParamAccess.item, 40);
            pManager.AddIntegerParameter("N", "N", "N", GH_ParamAccess.item, 3);
            pManager.AddCurveParameter("C", "C", "C", GH_ParamAccess.list);

            pManager[2].Optional = true;
            pManager[3].Optional = true;
            pManager[5].Optional = true;
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
            pManager.AddBrepParameter("Patches", "Patches", "Patches", GH_ParamAccess.list);
            pManager.AddNumberParameter("PlnError", "PlnError", "PlnError", GH_ParamAccess.item);

        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            int N = 0;
            m_Curves = new List<Curve>();
            var ExtCurve = false;

            if (DA.GetDataList(5, m_Curves))
                ExtCurve = true;
            else
            {
                m_Curves = new List<Curve>();
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
                m_Curves.Add(Crv1);
                m_Curves.Add(Crv2);
                m_Curves.Add(Crv3);
                m_Curves.Add(Crv4);
                m_Curves.Add(Crv5);
            }


            List<double> Dvalues = new List<double>();
            int Side = 0;

            //DA.GetDataList(0, Dvalues);
            //DA.GetData(1,ref Side);

            double u = 0, v = 0;
            int u_scl=0,  v_scl=0 ;
            DA.GetData(0, ref u);
            DA.GetData(1, ref v);
            DA.GetData(2, ref u_scl);
            DA.GetData(3, ref v_scl);
            DA.GetData(4, ref N);
            if (u_scl == 0 || v_scl == 0) return;

            double Kappa = 0;

            double landa = ComputeSideBLFunction(Dvalues, Side);
            if (Side != 0)
                Kappa = ComputeCornerBLFunction(Dvalues, Side - 1, Side);
            else
                Kappa = ComputeCornerBLFunction(Dvalues, Dvalues.Count - 1, Side);

            double Nu = ComputeSPsideBLFunction(Dvalues, Side);

            ComputeDomainPolygon();
            ComputePatches(N);

            //m_SrfPt = Kato_Suv(u, v);
            //ComputeSrfPts(u_scl, v_scl);
            var num_u = u_scl;
            var num_v = v_scl;
            ComputeSurfacePoints2(num_u, num_v);

            DA.SetDataList(0, m_Curves);
            DA.SetDataList(1, m_DPolygonLines);
            DA.SetData(2, landa);
            DA.SetData(3, Kappa);
            DA.SetData(4, Nu);
            DA.SetData(5, m_SrfPt);
            DA.SetDataList(6, m_SrfPts);
            //DA.SetDataList(7, m_Patches);
            DA.SetData(8, m_PlnAverageError);
        }

        private void ComputePatches(int N)
        {
            Point3d CenterPt = new Point3d();

            for (int i = 0; i < m_DomainPolygon.Count; i++)
                CenterPt = CenterPt + m_DomainPolygon[i];
            CenterPt = CenterPt / (double)m_DomainPolygon.Count;

            List<List<Point3d>> SubPolygons = new List<List<Point3d>>();
            List<Point3d> SubPolygon;
            Point3d P, P1, P2;
            int j, k;
            for (int i = 0; i < m_DomainPolygon.Count; i++)
            {
                SubPolygon = new List<Point3d>();
                j = (i + 1) % m_DomainPolygon.Count;
                k = (i - 1 + m_DomainPolygon.Count) % m_DomainPolygon.Count;
                P = m_DomainPolygon[i];
                P1 = 0.5 * (m_DomainPolygon[i] + m_DomainPolygon[j]);
                P2 = 0.5 * (m_DomainPolygon[i] + m_DomainPolygon[k]);

                SubPolygon.Add(P);
                SubPolygon.Add(P1);
                SubPolygon.Add(CenterPt);
                SubPolygon.Add(P2);

                SubPolygons.Add(SubPolygon);
            }

            List<List<List<Point3d>>> UVsss = new List<List<List<Point3d>>>();
            List<List<Point3d>> UVss;
            for (int i = 0; i < SubPolygons.Count; i++)
            {
                UVss = ComputeUVs(SubPolygons[i], N);
                UVsss.Add(UVss);
            }

            List<List<Point3d>> PatchUVs = new List<List<Point3d>>();
            List<Point3d> PatchUV;
            for (int i = 0; i < UVsss.Count; i++)
            {
                for (j = 0; j < UVsss[i].Count - 1; j++)
                {
                    for (k = 0; k < UVsss[i][j].Count - 1; k++)
                    {
                        PatchUV = new List<Point3d>();
                        PatchUV.Add(new Point3d(UVsss[i][j][k]));
                        PatchUV.Add(new Point3d(UVsss[i][j][k + 1]));
                        PatchUV.Add(new Point3d(UVsss[i][j + 1][k + 1]));
                        PatchUV.Add(new Point3d(UVsss[i][j + 1][k]));
                        PatchUVs.Add(PatchUV);
                    }
                }
            }

            SetBPatches(PatchUVs);
        }

        private void SetBPatches(List<List<Point3d>> PatchUVs)
        {
            m_Patches = new List<Brep>();
            Line ln1, ln2, ln3, ln4;
            List<Curve> lines = new List<Curve>();
            Brep BPatch;
            Point3d Pt1, Pt2, Pt3, Pt4;
            m_PlanarityErrors = new List<double>();
            double error;

            for (int j = 0; j < PatchUVs.Count; j++)
            {
                Pt1 = Kato_Suv(PatchUVs[j][0].X, PatchUVs[j][0].Y);
                Pt2 = Kato_Suv(PatchUVs[j][1].X, PatchUVs[j][1].Y);
                Pt3 = Kato_Suv(PatchUVs[j][2].X, PatchUVs[j][2].Y);
                Pt4 = Kato_Suv(PatchUVs[j][3].X, PatchUVs[j][3].Y);

                ln1 = new Line(Pt1, Pt2);
                ln2 = new Line(Pt2, Pt3);
                ln3 = new Line(Pt3, Pt4);
                ln4 = new Line(Pt4, Pt1);

                lines.Clear();
                lines.Add(ln1.ToNurbsCurve());
                lines.Add(ln2.ToNurbsCurve());
                lines.Add(ln3.ToNurbsCurve());
                lines.Add(ln4.ToNurbsCurve());

                BPatch = Brep.CreateEdgeSurface(lines);
                m_Patches.Add(BPatch);

                error = ComputePlanarityError(Pt1, Pt2, Pt3, Pt4);
                m_PlanarityErrors.Add(error);
            }

            m_PlnAverageError = 0;
            System.IO.File.Delete("C://result//Transfinite//Errors.csv");
            for (int i = 0; i < m_PlanarityErrors.Count; i++)
            {
                m_PlnAverageError = m_PlnAverageError + m_PlanarityErrors[i];
                System.IO.File.AppendAllText("C://result//Transfinite//Errors.csv", m_PlanarityErrors[i].ToString() + "\n");
            }
            m_PlnAverageError = m_PlnAverageError / (double)m_PlanarityErrors.Count;
        }

        private double ComputePlanarityError(Point3d Pt1, Point3d Pt2, Point3d Pt3, Point3d Pt4)
        {
            double error = 0;

            Plane pln1 = new Plane(Pt1, Pt2, Pt3);
            Plane pln2 = new Plane(Pt2, Pt3, Pt4);
            Plane pln3 = new Plane(Pt3, Pt4, Pt1);
            Plane pln4 = new Plane(Pt4, Pt1, Pt2);

            error = error + pln1.DistanceTo(Pt4);
            error = error + pln2.DistanceTo(Pt1);
            error = error + pln3.DistanceTo(Pt2);
            error = error + pln4.DistanceTo(Pt3);
            error = error / 4.0;

            return error;
        }

        private List<List<Point3d>> ComputeUVs(List<Point3d> SubPolygon, int N)
        {
            List<List<Point3d>> UVss = new List<List<Point3d>>();
            List<Point3d> UVs;
            Point3d UV;

            for (double v = 0; v < 1.00001; v = v + 1.0 / (double)N)
            {
                UVs = new List<Point3d>();
                for (double u = 0; u < 1.00001; u = u + 1.0 / (double)N)
                {
                    UV = new Point3d();
                    UV = (1 - u) * (1 - v) * SubPolygon[0] + u * (1 - v) * SubPolygon[1] + u * v * SubPolygon[2] + (1 - u) * v * SubPolygon[3];
                    UVs.Add(UV);
                }
                UVss.Add(UVs);
            }

            return UVss;
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
            m_uv = new List<Point3d>();
            var PlyPts = poly.ToNurbsCurve().Points;

            var Ulist = new List<double>();
            var Vlist = new List<double>();


            for (double i = m_minU; i <= m_maxU; i += inc_U)
            { Ulist.Add(i); }

            for (double j = m_minV; j <= m_maxV; j += inc_V)
            { Vlist.Add(j); }

            for (int i = 0; i < PlyPts.Count; i++)
            {

                if (!Ulist.Contains(PlyPts[i].X))
                    Ulist.Add(PlyPts[i].X);
                if (!Vlist.Contains(PlyPts[i].Y))
                    Vlist.Add(PlyPts[i].Y);
            }

            Ulist.Sort();
            Vlist.Sort();
            Ulist.RemoveAt(0);
            Ulist.RemoveAt(Ulist.Count - 1);
            Vlist.RemoveAt(0);
            Vlist.RemoveAt(Vlist.Count - 1);


            foreach (var UPt in Ulist)
            {
                var counterx = 0;
                var Ln = new Line(new Point3d(UPt, m_minV, 0), new Point3d(UPt, m_maxV, 0));
                var sr = Intersection.CurveLine(poly, Ln, 0.0001, 0.001);
                foreach (var VPt in Vlist)
                {
                    var state = poly.Contains(new Point3d(UPt, VPt, 0), Plane.WorldXY, 0.00001);
                    counterx++;
                    if (state != PointContainment.Outside)
                    {
                        //if (m_DomainPolygon.Contains(new Point3d(i, j, 0)))
                        //   m_SrfPts.Add(m_Curves[m_DomainPolygon.IndexOf(new Point3d(i, j, 0))].PointAtStart);
                        //else
                        m_SrfPts.Add(Kato_Suv(UPt, VPt));
                        m_uv.Add(new Point3d(UPt, VPt, 0));
                    }
                    else
                    {
                        var new_pt = new Point3d(UPt, VPt, 0);
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
                        m_uv.Add(new Point3d(edge_point[idx].X, edge_point[idx].Y, 0));
                    }

                }
                intlist.Add(counterx);
            }
        }

        private void ComputeSurfacePoints2(double u_div, double v_div)
        {
            m_SrfPts = new List<Point3d>();

            var polygon = new List<Point3d>(m_DomainPolygon);
            polygon.Add(new Point3d(polygon[0]));

            List<Line> PolyLine = new List<Line>();

            for (int i = 0; i < polygon.Count - 1; i++)
                PolyLine.Add(new Line(polygon[i], polygon[i + 1]));

            var u_inc = 1 / (u_div);
            v_div = (v_div / 2);

            var midLine = new Line(PolyLine[0].PointAtLength(PolyLine[0].Length / 2), PolyLine[2].To);
            m_DPolygonLines.Add(midLine);


            var DividingLines = new List<Line>();
            for (int i = 0; i <= u_div; i++)
            {
                DividingLines.Add(new Line(PolyLine[4].PointAtLength(PolyLine[4].Length - (double)i / u_div * PolyLine[4].Length), midLine.PointAtLength((double)i / u_div * midLine.Length)));
                DividingLines.Add(new Line(midLine.PointAtLength((double)i / u_div * midLine.Length), PolyLine[1].PointAtLength((double)i / u_div * PolyLine[1].Length)));
            }
            m_DPolygonLines.AddRange(DividingLines);


            for (int i = 0; i < DividingLines.Count / 2; i++)
            {
                for (int j = 0; j < v_div; j++)
                {
                    var pt = DividingLines[i * 2].PointAtLength((double)j / v_div * DividingLines[i * 2].Length);
                    m_SrfPts.Add(Kato_Suv(pt.X, pt.Y));
                    //m_SrfPts.Add(pt);

                }
                for (int j = 0; j <= v_div; j++)
                {
                    var pt = DividingLines[i * 2 + 1].PointAtLength((double)j / v_div * DividingLines[i * 2 + 1].Length);
                    m_SrfPts.Add(Kato_Suv(pt.X, pt.Y));
                    //m_SrfPts.Add(pt);

                }
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