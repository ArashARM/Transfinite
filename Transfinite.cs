using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using Rhino.Geometry.Intersect;
using System.Linq;
using System.Diagnostics;
using static GrasshopperProjects.Harmonic;//for harmonic
using Rhino.Commands;
using System.IO;
using Accord;
using Grasshopper;
using Grasshopper.Kernel.Data;
//using Accord.Math;
//using Accord.Math.Geometry;

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
        private static List<Curve> m_CurvesSerhat;
        private static List<Vector3d> m_Ts;
        private static List<List<Vector3d>> m_Tss;
        private static List<List<Vector3d>> m_TssSerhat;
        private static List<Curve> m_IsoCurves;
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
        private List<List<Curve>> IsolinesList = new List<List<Curve>>();

        public static List<Point3d> m_uvSerhat;

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
            pManager.AddGenericParameter("landa", "landa", "landa", GH_ParamAccess.item);
            pManager.AddGenericParameter("kappa", "kappa", "kappa", GH_ParamAccess.item);
            pManager.AddGenericParameter("nu", "nu", "nu", GH_ParamAccess.item);
            pManager.AddGenericParameter("kato", "kato", "kato", GH_ParamAccess.list);
            pManager.AddGenericParameter("kato_pointlist", "kato_pointlist", "kato_pointlist", GH_ParamAccess.list);
            pManager.AddBrepParameter("Patches", "Patches", "Patches", GH_ParamAccess.list);
            pManager.AddNumberParameter("PlnError", "PlnError", "PlnError", GH_ParamAccess.item);
            pManager.AddGenericParameter("Maplist", "Maplist", "Map string list", GH_ParamAccess.tree);

        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            int N = 0;
            m_Curves = new List<Curve>();
            m_CurvesSerhat = new List<Curve>();
            List<List<Point3d>> curvepoints = new List<List<Point3d>>();
            var ExtCurve = false;

            if (DA.GetDataList(5, m_Curves))
                ExtCurve = true;
            else
            {
                m_Curves = new List<Curve>();
                /*List<Point3d> Pts1 = new List<Point3d>() { new Point3d(0, 0, 0), new Point3d(3, 0, 1), new Point3d(6, 0, -1), new Point3d(10, -1, 0) };
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
                m_Curves.Add(Crv5);*/
                /* List<Point3d> Pts1 = new List<Point3d>() { new Point3d(0, 0, 0), new Point3d(3, 0, 1), new Point3d(6, 0, -1), new Point3d(10, 0, 0) };
                 List<Point3d> Pts2 = new List<Point3d>() { new Point3d(10, 0, 0), new Point3d(11, 5, -1), new Point3d(10, 8, 1), new Point3d(10, 10, 0) };
                 List<Point3d> Pts3 = new List<Point3d>() { new Point3d(10, 10, 0), new Point3d(8, 9, 1), new Point3d(6.5, 9, 0), new Point3d(5, 10, 0) };
                 List<Point3d> Pts4 = new List<Point3d>() { new Point3d(5, 10, 0), new Point3d(5, 14, 0), new Point3d(5, 18, 0), new Point3d(5, 20, 0) };
                 List<Point3d> Pts5 = new List<Point3d>() { new Point3d(5, 20, 0), new Point3d(4, 20, 1), new Point3d(2, 21, 0), new Point3d(0, 20, 0) };
                 List<Point3d> Pts6 = new List<Point3d>() { new Point3d(0, 20, 0), new Point3d(0, 14, 1), new Point3d(0, 8, -1), new Point3d(0, 0, 0) };*/
                /*List<Point3d> Pts1 = new List<Point3d>() { new Point3d(0, 0, 0), new Point3d(3, 0, 0), new Point3d(6, 0, 0), new Point3d(10, 0, 0) };
                List<Point3d> Pts2 = new List<Point3d>() { new Point3d(10, 0, 0), new Point3d(11, 5, 0), new Point3d(10, 8, 0), new Point3d(10, 10, 0) };
                List<Point3d> Pts3 = new List<Point3d>() { new Point3d(10, 10, 0), new Point3d(8, 9, 0), new Point3d(6.5, 9, 0), new Point3d(5, 10, 0) };
                List<Point3d> Pts4 = new List<Point3d>() { new Point3d(5, 10, 0), new Point3d(5, 14, 0), new Point3d(5, 18, 0), new Point3d(5, 20, 0) };
                List<Point3d> Pts5 = new List<Point3d>() { new Point3d(5, 20, 0), new Point3d(4, 20, 0), new Point3d(2, 21, 0), new Point3d(0, 20, 0) };
                List<Point3d> Pts6 = new List<Point3d>() { new Point3d(0, 20, 0), new Point3d(0, 14, 0), new Point3d(0, 8, 0), new Point3d(0, 0, 0) };*/
                List<Point3d> Pts1 = new List<Point3d>() { new Point3d(0, 0, 0), new Point3d(3, 0, 0), new Point3d(6, 0, 0), new Point3d(10, 0, 0) };
                List<Point3d> Pts2 = new List<Point3d>() { new Point3d(10, 0, 0), new Point3d(10, 5, 0), new Point3d(10, 8, 0), new Point3d(10, 10, 0) };
                List<Point3d> Pts3 = new List<Point3d>() { new Point3d(10, 10, 0), new Point3d(8, 10, 0), new Point3d(6.5, 10, 0), new Point3d(5, 10, 0) };
                List<Point3d> Pts4 = new List<Point3d>() { new Point3d(5, 10, 0), new Point3d(5, 14, 0), new Point3d(5, 18, 0), new Point3d(5, 20, 0) };
                List<Point3d> Pts5 = new List<Point3d>() { new Point3d(5, 20, 0), new Point3d(4, 20, 0), new Point3d(2, 20, 0), new Point3d(0, 20, 0) };
                List<Point3d> Pts6 = new List<Point3d>() { new Point3d(0, 20, 0), new Point3d(0, 14, 0), new Point3d(0, 8, 0), new Point3d(0, 0, 0) };

                /*List<Point3d> Pts1 = new List<Point3d>() { new Point3d(0, 0, 0), new Point3d(3, 0, 1), new Point3d(6, 0, -1), new Point3d(10, 0, 0) };
                List<Point3d> Pts2 = new List<Point3d>() { new Point3d(10, 0, 0), new Point3d(11, 5, -1), new Point3d(10, 8, 1), new Point3d(10, 10, 0) };
                List<Point3d> Pts3 = new List<Point3d>() { new Point3d(10, 10, 0), new Point3d(9, 12, 1), new Point3d(8.5, 15, 0), new Point3d(8, 18, 0) };
                List<Point3d> Pts4 = new List<Point3d>() { new Point3d(8, 18, 0), new Point3d(7, 19, 0), new Point3d(6, 19.5, 0), new Point3d(5, 20, 0) };
                List<Point3d> Pts5 = new List<Point3d>() { new Point3d(5, 20, 0), new Point3d(4, 20, 1), new Point3d(2, 21, 0), new Point3d(0, 20, 0) };
                List<Point3d> Pts6 = new List<Point3d>() { new Point3d(0, 20, 0), new Point3d(0, 14, 1), new Point3d(0, 8, -1), new Point3d(0, 0, 0) };*/

                Curve Crv1 = Curve.CreateInterpolatedCurve(Pts1, 3);
                Curve Crv2 = Curve.CreateInterpolatedCurve(Pts2, 3);
                Curve Crv3 = Curve.CreateInterpolatedCurve(Pts3, 3);
                Curve Crv4 = Curve.CreateInterpolatedCurve(Pts4, 3);
                for (int add = 1; add < Pts4.Count; add++)
                {
                    Pts3.Add(Pts4[add]);
                }
                Curve Crv34 = Curve.CreateInterpolatedCurve(Pts3, 1);
                Curve Crv5 = Curve.CreateInterpolatedCurve(Pts5, 3);
                Curve Crv6 = Curve.CreateInterpolatedCurve(Pts6, 3);
                m_Curves.Add(Crv1);
                m_Curves.Add(Crv2);
                m_Curves.Add(Crv3);
                m_Curves.Add(Crv4);
                m_Curves.Add(Crv5);
                m_Curves.Add(Crv6);

                curvepoints.Add(Pts1);
                curvepoints.Add(Pts2);
                curvepoints.Add(Pts3);
                curvepoints.Add(Pts5);
                curvepoints.Add(Pts6);

                m_CurvesSerhat.Add(Crv1);
                m_CurvesSerhat.Add(Crv2);
                m_CurvesSerhat.Add(Crv34);
                m_CurvesSerhat.Add(Crv5);
                m_CurvesSerhat.Add(Crv6);

                m_Ts = new List<Vector3d>();
                m_Ts.Add(new Vector3d(1, 1, 1));
                m_Ts.Add(new Vector3d(-1, 1, 1));
                m_Ts.Add(new Vector3d(-1, -1, 1));
                m_Ts.Add(new Vector3d(-1, -1, 1));
                m_Ts.Add(new Vector3d(-1, -1, 1));
                m_Ts.Add(new Vector3d(1, -1, 1));

                m_Tss = new List<List<Vector3d>>();
                List<Vector3d> Ts = new List<Vector3d>();
                Ts.Add(new Vector3d(0, 1, 1));
                Ts.Add(new Vector3d(0, 1, 1));
                m_Tss.Add(Ts);

                Ts = new List<Vector3d>();
                Ts.Add(new Vector3d(-1, 0, 1));
                Ts.Add(new Vector3d(-1, 0, 1));
                m_Tss.Add(Ts);

                Ts = new List<Vector3d>();
                Ts.Add(new Vector3d(0, -1, 1));
                Ts.Add(new Vector3d(0, -1, 1));
                m_Tss.Add(Ts);

                Ts = new List<Vector3d>();
                Ts.Add(new Vector3d(-1, 0, 1));
                Ts.Add(new Vector3d(-1, 0, 1));
                m_Tss.Add(Ts);

                Ts = new List<Vector3d>();
                Ts.Add(new Vector3d(0, -1, 1));
                Ts.Add(new Vector3d(0, -1, 1));
                m_Tss.Add(Ts);

                Ts = new List<Vector3d>();
                Ts.Add(new Vector3d(1, 0, 1));
                Ts.Add(new Vector3d(1, 0, 1));
                m_Tss.Add(Ts);

                ////

                m_TssSerhat = new List<List<Vector3d>>();
                List<Vector3d> TsSerhat = new List<Vector3d>();
                TsSerhat.Add(new Vector3d(1, 1, 0));
                TsSerhat.Add(new Vector3d(-1, 1, 0));
                m_TssSerhat.Add(TsSerhat);

                TsSerhat = new List<Vector3d>();
                TsSerhat.Add(new Vector3d(-1, 1, 1));
                TsSerhat.Add(new Vector3d(-1, -1, 1));
                m_TssSerhat.Add(TsSerhat);

                TsSerhat = new List<Vector3d>();
                TsSerhat.Add(new Vector3d(-1, -1, 1));
                TsSerhat.Add(new Vector3d(-1, -1, 1));
                m_TssSerhat.Add(TsSerhat);

                TsSerhat = new List<Vector3d>();
                TsSerhat.Add(new Vector3d(-1, -1, 1));
                TsSerhat.Add(new Vector3d(1, -1, 1));
                m_TssSerhat.Add(TsSerhat);

                //TsSerhat = new List<Vector3d>();
                //TsSerhat.Add(new Vector3d(0, -1, 1));
                //TsSerhat.Add(new Vector3d(0, -1, 1));
                //m_TssSerhat.Add(TsSerhat);

                TsSerhat = new List<Vector3d>();
                TsSerhat.Add(new Vector3d(1, -1, 1));
                TsSerhat.Add(new Vector3d(1, 1, 1));
                m_TssSerhat.Add(TsSerhat);
            }

            List<double> Dvalues = new List<double>();

            double u = 0, v = 0;
            int u_scl = 0, v_scl = 0;
            DA.GetData(0, ref u);
            DA.GetData(1, ref v);
            DA.GetData(2, ref u_scl);
            DA.GetData(3, ref v_scl);
            DA.GetData(4, ref N);
            if (u_scl == 0 || v_scl == 0) return;

            m_uvSerhat = new List<Point3d>();

            ComputeDomainPolygon();
            //HarmonicMapCreate_line(m_Curves.Count);
            //HarmonicMapCreate_curve(m_Curves.Count);
            //HarmonicMapCreate_curve1(m_Curves.Count);
            HarmonicMapCreate_curve(curvepoints, curvepoints.Count);

            //ComputeIsolines();



            Point3d Pt;
            Vector3d Vec;
            int k;
            double t;
            //for (int i = 0; i < m_Curves.Count; i++)
            //{
            //    k = (i + 1) % m_Curves.Count;
            //    for (int j = 0; j <= 100; j++)
            //    {
            //        t = m_Curves[i].Domain.Min + (m_Curves[i].Domain.Max - m_Curves[i].Domain.Min) * (double)j / 100.0;
            //        Pt = m_Curves[i].PointAt(t);
            //        //Vec = m_Ts[i] + (m_Ts[k] - m_Ts[i]) * (double)j/100.0;
            //        Vec = m_Tss[i][0] + (m_Tss[i][1] - m_Tss[i][0]) * (double)j / 100.0;


            //        m_DPolygonLines.Add(new Line(Pt, Pt + Vec));
            //    }
            //}

            ///*******
            # region ribbon visualization
            ////
            List<Line> ribbonlines = new List<Line>();
            for (int i = 0; i < m_CurvesSerhat.Count; i++)
            {
                k = (i + 1) % m_CurvesSerhat.Count;
                for (int j = 0; j <= 100; j++)
                {
                    t = m_CurvesSerhat[i].Domain.Min + (m_CurvesSerhat[i].Domain.Max - m_CurvesSerhat[i].Domain.Min) * (double)j / 100.0;
                    Pt = m_CurvesSerhat[i].PointAt(t);
                    //Vec = m_Ts[i] + (m_Ts[k] - m_Ts[i]) * (double)j/100.0;
                    Vec = m_TssSerhat[i][0] + (m_TssSerhat[i][1] - m_TssSerhat[i][0]) * (double)j / 100.0;

                    ribbonlines.Add(new Line(Pt, Pt + Vec));
                }
            }

            Point3d single_pt = new Point3d();
            Line single_vec = new Line();

            int curveidx = 2;
            t = m_CurvesSerhat[curveidx].Domain.Min + (m_CurvesSerhat[curveidx].Domain.Max - m_CurvesSerhat[curveidx].Domain.Min) * (u);
            single_pt = m_CurvesSerhat[curveidx].PointAt(t);
            //Vec = m_TssSerhat[curveidx][0] + (m_TssSerhat[curveidx][1] - m_TssSerhat[curveidx][0]) * (double)u / 100.0;
            Vec = m_TssSerhat[curveidx][0] + (u / (m_CurvesSerhat[curveidx].Domain.Max - m_CurvesSerhat[curveidx].Domain.Min)) * (m_TssSerhat[curveidx][1] - m_TssSerhat[curveidx][0]);
            single_vec = new Line(single_pt, single_pt + Vec);
            ////
            #endregion
            ///**********


            //ComputePatches(N);

            //m_SrfPt = Kato_Suv(u, v);
            //ComputeSrfPts(u_scl, v_scl);
            var num_u = u_scl;
            var num_v = v_scl;
            //ComputeSurfacePoints2(num_u, num_v);
            //ComputeSurfacePoints3(10);

            ComputeSurfacePoints4(10);

            //m_SrfPts.Clear();
            //m_SrfPts.Add(Kato_Suv(0, 0));
            //m_SrfPts.Add(Kato_Suv(1, 2));
            //m_SrfPts.Add(Kato_Suv(2, 4));
            //m_SrfPts.Add(Kato_Suv(3, 6));
            //m_SrfPts.Add(Kato_Suv(4, 8));
            //m_SrfPts.Add(Kato_Suv(4.99, 9.99));
            //m_SrfPts.Add(Kato_Suv(8, 10));
            //m_SrfPts.Add(Kato_Suv(6, 10));
            //m_SrfPts.Add(Kato_Suv(5, 10));
            //m_SrfPts.Add(Kato_Suv(5, 11));
            //m_SrfPts.Add(Kato_Suv(5, 12));
            //m_SrfPts.Add(Kato_Suv(5, 16));


            //ComputeDistance3(0.2, 0.6);
            //m_SrfPts = new List<Point3d>();
            //m_SrfPts.Add(Kato_Suv(0.5, 0.5));
            //m_SrfPts.Add(Kato_Suv(1, 0.5));
            //m_Curves = m_IsoCurves;


            ///************
            #region Harmonic map datatree contruction for visualization
            ///
            DataTree<string> stringtree = new DataTree<string>();
            for (int i = 0; i < HarmonicMapList_si.Count; i++)
            {
                stringtree.Add(writeroutput(HarmonicMapList_si[i]), new GH_Path(new int[] { 0, i }));
            }
            for (int i = 0; i < HarmonicMapList_di.Count; i++)
            {
                stringtree.Add(writeroutput(HarmonicMapList_di[i]), new GH_Path(new int[] { 1, i }));
            }
            ///
            #endregion
            ///************

            DA.SetDataList(0, m_Curves);
            DA.SetDataList(1, ribbonlines);
            DA.SetData(2, single_pt);
            DA.SetData(3, single_vec);
            DA.SetData(4, null);
            DA.SetDataList(5, m_uvSerhat);
            DA.SetDataList(6, m_SrfPts);
            DA.SetDataList(7, m_Patches);
            DA.SetData(8, m_PlnAverageError);
            DA.SetDataTree(9, stringtree);
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
            //System.IO.File.Delete("C://result//Transfinite//Errors.csv");
            //for (int i = 0; i < m_PlanarityErrors.Count; i++)
            //{
            //    m_PlnAverageError = m_PlnAverageError + m_PlanarityErrors[i];
            //    System.IO.File.AppendAllText("C://result//Transfinite//Errors.csv", m_PlanarityErrors[i].ToString() + "\n");
            //}
            //m_PlnAverageError = m_PlnAverageError / (double)m_PlanarityErrors.Count;
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

        private List<(double, double)> MVC(double u, double v)
        {
            List<double> w_i = new List<double>();
            List<double> s_i = new List<double>();

            var pt = new Point3d(u, v, 0);
            MVCweights(pt, out w_i);

            var d_i = w_i.Select(x => 1 - Math.Abs(x)).ToList();
            for (int i = 0; i < w_i.Count; i++)
            {
                if (w_i[i] == 0 && w_i[(i + 1) % w_i.Count] == 0)
                    s_i.Add(1);
                else
                    s_i.Add((double)d_i[i] / (double)(d_i[i] + d_i[(i + 1) % d_i.Count]));
            }
            for (int i = 0; i < d_i.Count; i++)
            {
                d_i[i] = (1 - Math.Abs(w_i[i] + (w_i[(i + 1) % w_i.Count] - w_i[i])) * s_i[i]);
            }
            List<(double, double)> Distance = new List<(double, double)>();

            for (int i = 0; i < d_i.Count; i++)
                Distance.Add((s_i[i], d_i[i]));

            return Distance;
        }

        private int MVCweights(Point3d queryCoord, out List<double> baryCoords)
        {
            var cageCoords = m_DomainPolygon;
            int nSize = cageCoords.Count;
            Debug.Assert(nSize != 0);

            double dx, dy;

            List<Vector3d> s = new List<Vector3d>(nSize);
            baryCoords = new List<double>(nSize);

            for (int i = 0; i < nSize; i++)
            {
                dx = cageCoords[i].X - queryCoord.X;
                dy = cageCoords[i].Y - queryCoord.Y;
                s.Add(new Vector3d(dx, dy, 0));
            }

            for (int i = 0; i < nSize; i++)
                baryCoords.Add(0.0);

            int ip, im;      // (i+1) and (i-1)
            double ri, rp, Ai, Di, dl, mu;  // Distance
            double eps = 10.0 * 2.22507e-308;

            // First check if any coordinates close to the cage point or
            // lie on the cage boundary. These are special cases.
            for (int i = 0; i < nSize; i++)
            {
                ip = (i + 1) % nSize;
                ri = Math.Sqrt(s[i].X * s[i].X + s[i].Y * s[i].Y);
                Ai = 0.5 * (s[i].X * s[ip].Y - s[ip].X * s[i].Y);
                Di = s[ip].X * s[i].X + s[ip].Y * s[i].Y;
                if (ri <= eps)
                {
                    baryCoords[i] = 1.0;
                    return 0;
                }
                else if (Math.Abs(Ai) <= 0 && Di < 0.0)
                {
                    dx = cageCoords[ip].X - cageCoords[i].X;
                    dy = cageCoords[ip].Y - cageCoords[i].Y;
                    dl = Math.Sqrt(dx * dx + dy * dy);
                    Debug.Assert(dl > eps);
                    dx = queryCoord.X - cageCoords[i].X;
                    dy = queryCoord.Y - cageCoords[i].Y;
                    mu = Math.Sqrt(dx * dx + dy * dy) / dl;
                    Debug.Assert(mu >= 0.0 && mu <= 1.0);
                    baryCoords[i] = 1.0 - mu;
                    baryCoords[ip] = mu;
                    return 0;
                }
            }

            // Page #12, from the paper
            List<double> tanalpha = new List<double>(nSize); // tan(alpha/2)
            for (int i = 0; i < nSize; i++)
            {
                ip = (i + 1) % nSize;
                im = (nSize - 1 + i) % nSize;
                ri = Math.Sqrt(s[i].X * s[i].X + s[i].Y * s[i].Y);
                rp = Math.Sqrt(s[ip].X * s[ip].X + s[ip].Y * s[ip].Y);
                Ai = 0.5 * (s[i].X * s[ip].Y - s[ip].X * s[i].Y);
                Di = s[ip].X * s[i].X + s[ip].Y * s[i].Y;
                tanalpha.Add((ri * rp - Di) / (2.0 * Ai));
            }

            // Equation #11, from the paper
            double wi, wsum = 0.0;
            for (int i = 0; i < nSize; i++)
            {
                im = (nSize - 1 + i) % nSize;
                ri = Math.Sqrt(s[i].X * s[i].X + s[i].Y * s[i].Y);
                wi = 2.0 * (tanalpha[i] + tanalpha[im]) / ri;
                wsum += wi;
                baryCoords[i] = wi;
            }

            if (Math.Abs(wsum) > 0.0)
            {
                for (int i = 0; i < nSize; i++)
                    baryCoords[i] /= wsum;
            }

            return 0;

        }

        private void ComputeSurfacePoints4(int N)
        {
            m_SrfPts = new List<Point3d>();

            //for (double u = 0; u < 1.00001; u = u + 1.0 / (double)N) //0, 1.00001 
            //{
            //    for (double v = 0; v < 2.00001; v = v + 1.0 / (double)N) //0, 2.00001
            //    {
            //        if (IsInDomainPolygon(u, v))
            //            m_SrfPts.Add(Kato_Suv(u, v)); //m_SrfPts.Add(new Point3d(u, v, 0));
            //    }
            //}

            //for (double u = -20; u < 20.00001; u = u + 1.0) //0, 1.00001 
            //{
            //    for (double v = -20; v < 20.00001; v = v + 1.0) //0, 2.00001
            //    {
            //        if (IsInDomainPolygon(u, v))
            //            m_SrfPts.Add(new Point3d(u, v, 0)); //m_SrfPts.Add(Kato_Suv(u, v)); 
            //    }
            //}

            List<List<Point3d>> Ptss = new List<List<Point3d>>();
            List<Point3d> Pts;
            double v = -1;
            for (v = 0; v <= 10.0001; v = v + 1)
            {
                Pts = new List<Point3d>();
                for (double u = 0; u <= 10.0001; u = u + 1)
                {
                    //m_SrfPts.Add(Kato_Suv(u, v)); //m_SrfPts.Add(new Point3d(u, v, 0));
                    //Pts.Add(Kato_Suv(u, v));

                    Point3d onecalculation = Kato_Suv(u, v);
                    m_SrfPts.Add(onecalculation); //m_SrfPts.Add(new Point3d(u, v, 0));
                    Pts.Add(onecalculation);
                }
                Ptss.Add(Pts);
            }

            for (v = 11; v <= 20.0001; v = v + 1)
            {
                Pts = new List<Point3d>();
                for (double u = 0; u <= 10.0001; u = u + 1)
                {
                    if (u > 5)
                    {
                        //m_SrfPts.Add(new Point3d(u, v, 0));
                        Pts.Add(new Point3d(u, v, 0));
                    }
                    else
                    {
                        Point3d onecalculation = Kato_Suv(u, v);
                        //m_SrfPts.Add(Kato_Suv(u, v)); //m_SrfPts.Add(new Point3d(u, v, 0));
                        //Pts.Add(Kato_Suv(u, v));
                        m_SrfPts.Add(onecalculation); //m_SrfPts.Add(new Point3d(u, v, 0));
                        Pts.Add(onecalculation);
                    }
                }
                Ptss.Add(Pts);
            }

        }

        private void ComputeSurfacePoints3(int N)
        {
            m_SrfPts = new List<Point3d>();
            Line Ln = new Line(m_DomainPolygon[0], m_DomainPolygon[1]);
            Ln.Extend(100, 100);
            List<List<Point3d>> Ptss = new List<List<Point3d>>();
            List<Point3d> Pts = new List<Point3d>();
            Vector3d OffsetVec, Vec1, Vec2;

            for (int i = 0; i <= N; i++)
                Pts.Add(m_DomainPolygon[0] + i * (m_DomainPolygon[1] - m_DomainPolygon[0]) / (double)N);
            Ptss.Add(Pts);

            Vec1 = m_DomainPolygon[2] - m_DomainPolygon[1];
            Vec2 = m_DomainPolygon[m_DomainPolygon.Count - 1] - m_DomainPolygon[0];
            OffsetVec = 0.5 * (Vec1 + Vec2);
            OffsetVec.Unitize();
            OffsetVec = Vector3d.Multiply(0.1, OffsetVec);

            Ln.From = Ln.From + OffsetVec;
            Ln.To = Ln.To + OffsetVec;

            List<Line> Lines = new List<Line>();
            Line Ln2;
            int j;
            for (int i = 1; i < m_DomainPolygon.Count; i++)
            {
                //if (i == 2 || i == 4)
                //continue;
                j = (i + 1) % m_DomainPolygon.Count;
                Ln2 = new Line(m_DomainPolygon[i], m_DomainPolygon[j]);
                Lines.Add(Ln2);
            }

            List<Point3d> Pts2;
            Pts = GetIntersection(Ln, Lines);
            while (Pts.Count != 0)
            {
                if (Pts.Count <= 2)
                {
                    Ln.From = Ln.From + OffsetVec;
                    Ln.To = Ln.To + OffsetVec;
                    Pts = GetIntersection(Ln, Lines);
                    break;
                }
                Pts2 = new List<Point3d>();
                for (int i = 0; i <= N; i++)
                    Pts2.Add(Pts[0] + i * (Pts[1] - Pts[0]) / (double)N);
                Ptss.Add(Pts2);

                Ln.From = Ln.From + OffsetVec;
                Ln.To = Ln.To + OffsetVec;
                Pts = GetIntersection(Ln, Lines);



                //break;
            }

            for (int i = 0; i < Ptss.Count; i++)
                for (int k = 0; k < Ptss[i].Count; k++)
                    m_SrfPts.Add(Kato_Suv(Ptss[i][k].X, Ptss[i][k].Y));
            /*for (int i = 0; i < Ptss.Count; i++)
                for (int k = 0; k < Ptss[i].Count; k++)
                    m_SrfPts.Add(Ptss[i][k]);*/
            /*m_SrfPts.Add(new Point3d(0, 0, 0));
            m_SrfPts.Add(new Point3d(0.5, 0.5, 0));
            m_SrfPts.Add(new Point3d(0.75, 0.5, 0));
            m_SrfPts.Add(new Point3d(1, 0.5, 0));
            m_SrfPts.Add(new Point3d(0.25, 1.5, 0));*/
            /*m_SrfPts.Add(Kato_Suv(0, 0));
            m_SrfPts.Add(Kato_Suv(0.5, 0.5));
            m_SrfPts.Add(Kato_Suv(0.75, 0.5));
            m_SrfPts.Add(Kato_Suv(1, 0.5));
            m_SrfPts.Add(Kato_Suv(0.25, 1.5));*/

            //m_SrfPts.Add(Kato_Suv(0.5, 0.8));
            //m_SrfPts.Add(Kato_Suv(0.5, 0.5));
            //m_SrfPts.Add(Kato_Suv(0.5, 0.1));
            //m_SrfPts.Add(Kato_Suv(0.5, 1.5));
            //m_SrfPts.Add(Kato_Suv(0.2, 1));
            //m_SrfPts.Add(Kato_Suv(0.4, 1));
        }

        List<Point3d> GetIntersection(Line Ln, List<Line> Lines)
        {
            List<Point3d> Pts = new List<Point3d>();
            double a, b;

            for (int i = 0; i < Lines.Count; i++)
            {
                Intersection.LineLine(Ln, Lines[i], out a, out b);

                if (a >= 0 && a <= 1 && b >= 0 && b <= 1)
                {
                    Pts.Add(Ln.PointAt(a));
                }
            }

            return Pts;
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
            var ptt = new List<Point3d>();


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

            /*Line Ln = new Line(new Point3d(u, v, 0), new Point3d(u, v, 0) + new Vector3d(1000, 0, 0));
            int j;
            List<Line> Lines = new List<Line>();
            Line Ln2;
            for (int i = 0; i < m_DomainPolygon.Count; i++)
            {
                j = (i + 1) % m_DomainPolygon.Count;
                Ln2 = new Line(m_DomainPolygon[i], m_DomainPolygon[j]);
                Lines.Add(Ln2);
            }

            List<Point3d> Pts = GetIntersection(Ln, Lines);
            if (Pts.Count != 0)
                return true;  */

            return false;
        }

        public List<HarmonicMap> HarmonicMapList_si = new List<HarmonicMap>();
        public List<HarmonicMap> HarmonicMapList_di = new List<HarmonicMap>();
        private void HarmonicMapCreate_line(int n)
        {
            HarmonicMapList_si = new List<HarmonicMap>();
            HarmonicMapList_di = new List<HarmonicMap>();
            List<(double, double)> output = new List<(double, double)>();
            int i_before = 0, i_after = 0, i_afterafter = 0, i_afterafterafter = 0;
            int j;
            double value = 1;
            for (int i = 0; i < n; i++)
            {
                i_before = IndexWrapper((i - 1), n);
                i_after = IndexWrapper((i + 1), n);
                i_afterafter = IndexWrapper((i + 2), n);
                i_afterafterafter = IndexWrapper((i + 3), n);

                ///***********
                ///harmonic map created for si -----------------------------------
                ///***********
                double[] min = { 0, 0 }, max = { 20, 20 }; // domain polygon size
                HarmonicMap map_si;
                int levels = 9;//???
                map_si = harmonic_create(min, max, levels);

                ///Assigned Value on domain points
                List<Point3d> DomainPolygonHarmonic = new List<Point3d>(m_DomainPolygon);
                for (int ii = 0; ii < DomainPolygonHarmonic.Count; ii++)
                {
                    if (ii == i_after) DomainPolygonHarmonic[ii] = (new Point3d(DomainPolygonHarmonic[ii].X, DomainPolygonHarmonic[ii].Y, value));
                    else if (ii == i_afterafter) DomainPolygonHarmonic[ii] = (new Point3d(DomainPolygonHarmonic[ii].X, DomainPolygonHarmonic[ii].Y, value));
                    //else if (ii == i_afterafterafter) DomainPolygonHarmonic[ii] = (new Point3d(DomainPolygonHarmonic[ii].X, DomainPolygonHarmonic[ii].Y, value));
                    else if (ii == i_before) DomainPolygonHarmonic[ii] = (new Point3d(DomainPolygonHarmonic[ii].X, DomainPolygonHarmonic[ii].Y, 0));
                    else if (ii == i) DomainPolygonHarmonic[ii] = (new Point3d(DomainPolygonHarmonic[ii].X, DomainPolygonHarmonic[ii].Y, 0));
                    else DomainPolygonHarmonic[ii] = (new Point3d(DomainPolygonHarmonic[ii].X, DomainPolygonHarmonic[ii].Y, value));
                }
                ////

                for (int ii = 0; ii < DomainPolygonHarmonic.Count; ++ii)
                {
                    int j0 = (ii + 1) % 6;
                    harmonic_add_line(map_si, DomainPolygonHarmonic[ii], DomainPolygonHarmonic[j0]);
                }
                harmonic_solve(map_si, 1.0e-5, false);
                ///****
                HarmonicMapList_si.Add(map_si);
                ////
                ////--------------------------------------------------------------

                ///***********
                ///harmonic map created for di -----------------------------------
                ///***********
                double[] min_di = { 0, 0 }, max_di = { 20, 20 }; // domain polygon size
                HarmonicMap map_di;
                int levels_di = 9;//???
                map_di = harmonic_create(min_di, max_di, levels_di);

                ///Assigned Value on domain points
                DomainPolygonHarmonic = new List<Point3d>(m_DomainPolygon);
                for (int ii = 0; ii < DomainPolygonHarmonic.Count; ii++)
                {
                    if (ii == i) DomainPolygonHarmonic[ii] = (new Point3d(DomainPolygonHarmonic[ii].X, DomainPolygonHarmonic[ii].Y, 0));
                    else if (ii == i_after) DomainPolygonHarmonic[ii] = (new Point3d(DomainPolygonHarmonic[ii].X, DomainPolygonHarmonic[ii].Y, 0));
                    else DomainPolygonHarmonic[ii] = (new Point3d(DomainPolygonHarmonic[ii].X, DomainPolygonHarmonic[ii].Y, value));
                }
                ////

                for (int ii = 0; ii < DomainPolygonHarmonic.Count; ++ii)
                {
                    int j0 = (ii + 1) % 6;
                    harmonic_add_line(map_di, DomainPolygonHarmonic[ii], DomainPolygonHarmonic[j0]);
                }
                harmonic_solve(map_di, 1.0e-5, false);
                ///*****
                HarmonicMapList_di.Add(map_di);
                ////
                ////-------------------------------------------------------------
            }

        }

        private void HarmonicMapCreate_curve(List<List<Point3d>> curvepointlist, int n)
        {
            HarmonicMapList_si = new List<HarmonicMap>();
            HarmonicMapList_di = new List<HarmonicMap>();
            List<(double, double)> output = new List<(double, double)>();
            int i_before, i_after;
            double value = 1;

            //// Domain curve creation without "value" or "z"
            List<Point3d[]> DomainPolygonHarmonic_main = new List<Point3d[]> { };
            for (int i = 0; i < n; i++)
            {
                DomainPolygonHarmonic_main.Add(curvepointlist[i].ToArray());
            }

            for (int i = 0; i < DomainPolygonHarmonic_main.Count; i++)
            {
                i_before = IndexWrapper((i - 1), DomainPolygonHarmonic_main.Count);
                int i_beforebefore = IndexWrapper((i - 2), DomainPolygonHarmonic_main.Count);
                i_after = IndexWrapper((i + 1), DomainPolygonHarmonic_main.Count);

                ///***********
                ///harmonic map created for si -----------------------------------
                ///***********
                double[] min = { 0, 0 }, max = { 20, 20 }; // domain polygon size
                HarmonicMap map_si;
                int levels = 9;//???
                map_si = harmonic_create(min, max, levels);

                ///Assigned Value on domain points, si
                List<Point3d[]> DomainPolygonHarmonic = new List<Point3d[]>();
                for (int ii = 0; ii < DomainPolygonHarmonic_main.Count; ii++)
                {
                    DomainPolygonHarmonic.Add(new Point3d[DomainPolygonHarmonic_main[ii].Length]);
                    if (ii == i_after)
                    {
                        for (int pt = 0; pt < DomainPolygonHarmonic[ii].Length; pt++)
                        {
                            DomainPolygonHarmonic[ii][pt] = new Point3d(DomainPolygonHarmonic_main[ii][pt].X, DomainPolygonHarmonic_main[ii][pt].Y, value);
                        }
                    }
                    else if (ii == i_before)
                    {
                        for (int pt = 0; pt < DomainPolygonHarmonic[ii].Length; pt++)
                        {
                            DomainPolygonHarmonic[ii][pt] = new Point3d(DomainPolygonHarmonic_main[ii][pt].X, DomainPolygonHarmonic_main[ii][pt].Y, 0);
                        }
                    }
                    else if (ii == i)
                    {
                        double stepsize = value / (DomainPolygonHarmonic[ii].Length - 1);
                        for (int pt = 0; pt < DomainPolygonHarmonic[ii].Length; pt++)
                        {
                            DomainPolygonHarmonic[ii][pt] = new Point3d(DomainPolygonHarmonic_main[ii][pt].X, DomainPolygonHarmonic_main[ii][pt].Y, stepsize * pt);
                        }
                    }
                    else if (ii == i_beforebefore)
                    {
                        double stepsize = value / (DomainPolygonHarmonic[ii].Length - 1);
                        for (int pt = 0; pt < DomainPolygonHarmonic[ii].Length; pt++)
                        {
                            DomainPolygonHarmonic[ii][pt] = new Point3d(DomainPolygonHarmonic_main[ii][pt].X, DomainPolygonHarmonic_main[ii][pt].Y, stepsize * ((DomainPolygonHarmonic[ii].Length - 1) - pt));
                        }
                    }
                    else
                    {
                        for (int pt = 0; pt < DomainPolygonHarmonic[ii].Length; pt++)
                        {
                            DomainPolygonHarmonic[ii][pt] = new Point3d(DomainPolygonHarmonic_main[ii][pt].X, DomainPolygonHarmonic_main[ii][pt].Y, value);
                        }
                    }
                }
                ////

                for (int ii = 0; ii < DomainPolygonHarmonic.Count; ++ii)
                {
                    ///
                    ///middle point repeat
                    if (ii == 2)
                    {
                        List<Point3d> idle = DomainPolygonHarmonic[ii].ToList();
                        for (int xx = 0; xx < idle.Count; xx++)
                        {
                            if (xx == 3)
                            {
                                Point3d asd = new Point3d(idle[xx]);
                                idle.Insert(xx, asd);
                                idle.Insert(xx, asd);
                                idle.Insert(xx, asd);
                                idle.Insert(xx, asd);
                                //idle.Insert(xx, asd);
                                break;
                            }
                        }
                        DomainPolygonHarmonic[ii] = idle.ToArray();
                    }
                    ////
                    ///
                    harmonic_add_curve(map_si, DomainPolygonHarmonic[ii], DomainPolygonHarmonic[ii].Length);
                }
                harmonic_solve(map_si, 1.0e-5, false);
                ///****
                HarmonicMapList_si.Add(map_si);
                ////
                ////--------------------------------------------------------------

                ///***********
                ///harmonic map created for di -----------------------------------
                ///***********
                double[] min_di = { 0, 0 }, max_di = { 20, 20 }; // domain polygon size
                HarmonicMap map_di;
                int levels_di = 9;//???
                map_di = harmonic_create(min_di, max_di, levels_di);

                ///Assigned Value on domain points, di
                DomainPolygonHarmonic = new List<Point3d[]>();
                for (int ii = 0; ii < DomainPolygonHarmonic_main.Count; ii++)
                {
                    DomainPolygonHarmonic.Add(new Point3d[DomainPolygonHarmonic_main[ii].Length]);
                    if (ii == i_after)
                    {
                        double stepsize = value / (DomainPolygonHarmonic[ii].Length - 1);
                        for (int pt = 0; pt < DomainPolygonHarmonic[ii].Length; pt++)
                        {
                            DomainPolygonHarmonic[ii][pt] = new Point3d(DomainPolygonHarmonic_main[ii][pt].X, DomainPolygonHarmonic_main[ii][pt].Y, stepsize * pt);
                        }
                    }
                    else if (ii == i_before)
                    {
                        double stepsize = value / (DomainPolygonHarmonic[ii].Length - 1);
                        for (int pt = 0; pt < DomainPolygonHarmonic[ii].Length; pt++)
                        {
                            DomainPolygonHarmonic[ii][pt] = new Point3d(DomainPolygonHarmonic_main[ii][pt].X, DomainPolygonHarmonic_main[ii][pt].Y, stepsize * ((DomainPolygonHarmonic[ii].Length - 1) - pt));
                        }
                    }
                    else if (ii == i)
                    {
                        for (int pt = 0; pt < DomainPolygonHarmonic[ii].Length; pt++)
                        {
                            DomainPolygonHarmonic[ii][pt] = new Point3d(DomainPolygonHarmonic_main[ii][pt].X, DomainPolygonHarmonic_main[ii][pt].Y, 0);
                        }
                    }
                    else
                    {
                        for (int pt = 0; pt < DomainPolygonHarmonic[ii].Length; pt++)
                        {
                            DomainPolygonHarmonic[ii][pt] = new Point3d(DomainPolygonHarmonic_main[ii][pt].X, DomainPolygonHarmonic_main[ii][pt].Y, value);
                        }
                    }
                }
                ////

                for (int ii = 0; ii < DomainPolygonHarmonic.Count; ++ii)
                {
                    ///
                    ///middle point repeat
                    if (ii == 2)
                    {
                        List<Point3d> idle = DomainPolygonHarmonic[ii].ToList();
                        for (int xx = 0; xx < idle.Count; xx++)
                        {
                            if (xx == 3)
                            {
                                Point3d asd = new Point3d(idle[xx]);
                                idle.Insert(xx, asd);
                                idle.Insert(xx, asd);
                                idle.Insert(xx, asd);
                                idle.Insert(xx, asd);
                                idle.Insert(xx, asd);
                                break;
                            }
                        }
                        DomainPolygonHarmonic[ii] = idle.ToArray();
                    }
                    ////
                    ///
                    harmonic_add_curve(map_di, DomainPolygonHarmonic[ii], DomainPolygonHarmonic[ii].Length);
                }
                harmonic_solve(map_di, 1.0e-5, false);
                ///*****
                HarmonicMapList_di.Add(map_di);
                ////
                ////-------------------------------------------------------------
            }

        }

        private List<(double, double)> HarmonicMapCall(double u, double v, int n)
        {
            List<(double, double)> output = new List<(double, double)>();
            int j;
            Point3d point = new Point3d(u, v, 0);
            for (int i = 0; i < n; i++)
            {
                ///
                bool success_si = harmonic_eval(HarmonicMapList_si[i], point, out double result_si);//harmonic calculation
                ///
                bool success_di = harmonic_eval(HarmonicMapList_di[i], point, out double result_di);//harmonic calculation
                ///

                if (success_si == false) AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Harmonic function si!!");
                if (success_di == false) AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Harmonic function di!!");

                output.Add((result_si, result_di));//si_di
            }
            return output;
        }

        private Point3d Kato_Suv(double u, double v)
        {
            //List<(double, double)> si_di = ComputeRadialDistanceFunctionforKato(u, v, m_Curves.Count);
            //List<(double, double)> si_di = ComputeDistance2(u, v);
            //List<(double, double)> si_di = ComputeDistance3(u, v);
            //List<(double, double)> si_di = MVC(u, v);
            //List<(double, double)> si_di = HarmonicMapCall(u, v, m_Curves.Count);
            List<(double, double)> si_di = HarmonicMapCall(u, v, m_CurvesSerhat.Count);//curve part
            m_uvSerhat.Add(new Point3d(u, v, 0));
            //int j;
            //Line Ln;
            //Point3d Pt;
            //double s, d;
            /* for (int i = 0; i < m_DomainPolygon.Count; i++)
             {
                 j = (i+1)% m_DomainPolygon.Count;
                 Ln = new Line(m_DomainPolygon[i], m_DomainPolygon[j]);
                 Pt = Ln.ClosestPoint(new Point3d(u, v, 0), true);
                 s = Pt.DistanceTo(Ln.From) / Ln.Length;
                 d = Pt.DistanceTo(new Point3d(u, v, 0));

                 si_di.Add((s, d));
             }*/
            /*List<(double, double)> si_di = new List<(double, double)>();
             * si_di.Add((0.5, 0.8));
            si_di.Add((0.8, 0.5));
            si_di.Add((1, 0.2));
            si_di.Add((0.5, 1000.0));
            si_di.Add((0.5, 1000));
            si_di.Add((0.6, 1000));*/

            double s;
            List<double> d_i = new List<double>();
            for (int i = 0; i < si_di.Count; i++)
            {
                d_i.Add(si_di[i].Item2);
            }

            Vector3d T = new Vector3d();
            Point3d r_sum = new Point3d();
            List<double> Value = ComputeSPsideBLFunction1(d_i);
            int j;
            for (int i = 0; i < m_CurvesSerhat.Count; i++)
            {
                //double domainlengt = (m_DomainPolygon[IndexWrapper(i, m_Curves.Count)] - m_DomainPolygon[IndexWrapper(i + 1, m_Curves.Count)]).Length;
                double curveparameter = (m_CurvesSerhat[i].Domain.Length * si_di[i].Item1) / m_CurvesSerhat[i].GetLength();

                //Vector3d crossproduct = Vector3d.CrossProduct(m_Curves[i].TangentAt(curveparameter), m_Curves[i].CurvatureAt(curveparameter));
                //Point3d r = m_Curves[i].PointAt(curveparameter) + (si_di[i].Item2 * crossproduct);
                //s = m_Curves[i].Domain.Min + si_di[i].Item1 * (m_Curves[i].Domain.Max - m_Curves[i].Domain.Min);
                //Vector3d crossproduct = Vector3d.CrossProduct(m_Curves[i].TangentAt(s), m_Curves[i].CurvatureAt(s));
                //Vector3d crossproduct = m_Curves[i].CurvatureAt(s);
                //Point3d r = m_Curves[i].PointAt(s) + (si_di[i].Item2 * crossproduct);

                //Point3d r = m_Curves[i].PointAt(curveparameter) + (si_di[i].Item2 * m_Curves[i].TangentAt(curveparameter));
                //   s = m_Curves[i].Domain.Min + si_di[i].Item1 * (m_Curves[i].Domain.Max - m_Curves[i].Domain.Min);
                //   Point3d r = m_Curves[i].PointAt(s) + (si_di[i].Item2 * m_Curves[i].TangentAt(s));
                //double blendingfunctionvalue = ComputeSPsideBLFunction(d_i, i);

                s = m_CurvesSerhat[i].Domain.Min + si_di[i].Item1 * (m_CurvesSerhat[i].Domain.Max - m_CurvesSerhat[i].Domain.Min);
                j = (i + 1) % m_CurvesSerhat.Count;
                //T = m_Ts[i] + si_di[i].Item1 * (m_Ts[j] - m_Ts[i]);
                //if (m_TssSerhat[i].Count > 2)
                //{
                //    ;
                //}
                //else T = m_TssSerhat[i][0] + si_di[i].Item1 * (m_TssSerhat[i][1] - m_TssSerhat[i][0]);

                T = m_TssSerhat[i][0] + (si_di[i].Item1/ (m_CurvesSerhat[i].Domain.Max - m_CurvesSerhat[i].Domain.Min)) * (m_TssSerhat[i][1] - m_TssSerhat[i][0]);


                Point3d r = m_CurvesSerhat[i].PointAt(s) + (si_di[i].Item2 * T);

                r_sum += r * Value[i];
            }

            return r_sum;
        }

        private void ComputeIsolines()
        {
            IsolinesList.Clear();
            var Pnum = m_DPolygonLines.Count;
            for (int i = 0; i < Pnum; i++)
            {
                var lft = new Line();

                if (i == 0)
                    lft = m_DPolygonLines[Pnum - 1];
                else if (i == 2 || i == 5)
                    lft = m_DPolygonLines[(i + 1) % Pnum];
                else
                    lft = m_DPolygonLines[i - 1];


                var Left = lft.ToNurbsCurve();
                var Rght = new List<NurbsCurve>();
                if (i == 2 || i == 5)
                {
                    Rght.Add(m_DPolygonLines[(i - 1) % Pnum].ToNurbsCurve());
                    Rght.Add(m_DPolygonLines[(i - 2) % Pnum].ToNurbsCurve());
                    if (i == 2)
                        Rght.Add(m_DPolygonLines[Pnum - 1].ToNurbsCurve());
                    else
                        Rght.Add(m_DPolygonLines[(i - 3) % Pnum].ToNurbsCurve());
                }
                else
                {
                    Rght.Add(m_DPolygonLines[(i + 1) % Pnum].ToNurbsCurve());
                    Rght.Add(m_DPolygonLines[(i + 2) % Pnum].ToNurbsCurve());
                    Rght.Add(m_DPolygonLines[(i + 3) % Pnum].ToNurbsCurve());

                }


                var Right = Curve.JoinCurves(Rght);
                if (i == 0)
                {
                    //var rev2 = Left.Reverse();
                    var rev = Right[0].Reverse();
                }
                else if (i == 1)
                {
                    var rev2 = Left.Reverse();
                    //var rev = Right[0].Reverse();
                }
                else if (i == 2)
                {
                    var rev2 = Left.Reverse();
                    //var rev = Right[0].Reverse();
                }
                else if (i == 3)
                {
                    //var rev2 = Left.Reverse();
                    var rev = Right[0].Reverse();
                }
                else if (i == 4)
                {
                    var rev2 = Left.Reverse();
                    //var rev = Right[0].Reverse();
                }
                else if (i == 5)
                {
                    var rev2 = Left.Reverse();
                    //var rev = Right[0].Reverse();
                }

                IsolinesList.Add(Curve.CreateTweenCurves(Left, Right[0], 30, 0.0001).ToList());
                IsolinesList[i].Insert(0, Left);
                IsolinesList[i].Add(Right[0]);
                if (i == 2 || i == 5)
                    IsolinesList[i].Reverse();
                if (i == 0 || i == 2 || i == 3 || i == 5)
                    IsolinesList[i].ForEach(x => x.Reverse());
            }
            m_IsoCurves = IsolinesList[5];
        }

        private List<(double, double)> ComputeDistance(double u, double v)
        {
            List<(double, double)> Distance = new List<(double, double)>();
            var pt = new Point3d(u, v, 0);
            var tempdist = 0.0;
            var jlist = new List<int>();
            var slist = new List<double>();

            for (int i = 0; i < IsolinesList.Count; i++)
            {
                var dist = 9e9;
                int minj = 0;
                double t = 0;
                double treal = 0;
                for (int j = 0; j < IsolinesList[i].Count; j++)
                {
                    IsolinesList[i][j].ClosestPoint(pt, out t);
                    tempdist = pt.DistanceTo(IsolinesList[i][j].PointAt(t));
                    if (dist > tempdist)
                    {
                        dist = tempdist;
                        minj = j;
                        treal = t;
                    }

                }
                if (i == 2 || i == 5)
                    //Distance.Add((IsolinesList[i][minj].GetLength(new Interval(0,treal)), (minj / (IsolinesList[i].Count - 1))));
                    Distance.Add((((double)minj / (IsolinesList[i].Count - 1)), (IsolinesList[i][minj].GetLength(new Interval(treal, 1)))));
                //Distance.Add((((double)minj / (IsolinesList[i].Count - 1)), (IsolinesList[i][minj].GetLength(new Interval(0, treal)))));

                else
                    //Distance.Add((IsolinesList[i][minj].GetLength(new Interval(treal, 1)), (minj / (IsolinesList[i].Count - 1))));
                    Distance.Add((((double)minj / (IsolinesList[i].Count - 1)), (IsolinesList[i][minj].GetLength(new Interval(0, treal)))));
                //Distance.Add((((double)minj / (IsolinesList[i].Count - 1)), (IsolinesList[i][minj].GetLength(new Interval(treal, 1)))));

            }

            return Distance;
        }

        private List<(double, double)> ComputeDistance2(double u, double v)
        {
            List<(double, double)> Distance = new List<(double, double)>();
            double s = -1, d = -1, MinDt, Dt;
            Point3d Pt = new Point3d(u, v, 0), ClosestPt;
            double t;
            Interval subdomain;
            Point3d PtS, PtE;

            for (int i = 0; i < IsolinesList.Count; i++)
            {
                MinDt = 9E9;
                for (int j = 0; j < IsolinesList[i].Count; j++)
                {
                    PtS = IsolinesList[i][j].PointAtStart;
                    PtE = IsolinesList[i][j].PointAtEnd;
                    IsolinesList[i][j].ClosestPoint(Pt, out t);
                    ClosestPt = IsolinesList[i][j].PointAt(t);
                    Dt = Pt.DistanceTo(ClosestPt);

                    if (Dt < MinDt)
                    {
                        MinDt = Dt;

                        subdomain = new Interval(IsolinesList[i][j].Domain.Min, t);
                        d = IsolinesList[i][j].GetLength(subdomain);
                        s = (double)j / (double)(IsolinesList[i].Count - 1);
                    }
                }
                Distance.Add((s, d));
            }

            return Distance;
        }

        private List<(double, double)> ComputeDistance3(double u, double v)
        {
            List<(double, double)> Distance = new List<(double, double)>();
            int j, k;
            Vector3d Vec = new Vector3d(), Vec1 = new Vector3d(), Vec2 = new Vector3d(), Vec3 = new Vector3d();
            Point3d Pt = new Point3d(u, v, 0);
            List<double> ds = new List<double>(), ss = new List<double>();

            for (int i = 0; i < m_DomainPolygon.Count; i++)
            {
                j = (i + 1) % m_DomainPolygon.Count;
                Vec1 = Pt - m_DomainPolygon[i];
                Vec2 = Pt - m_DomainPolygon[j];
                Vec3 = m_DomainPolygon[j] - m_DomainPolygon[i];

                //Vec = Vec1 + Vec2 - Vec3;
                //ds.Add(Vec.Length);
                ds.Add(Vec1.Length + Vec2.Length - Vec3.Length);
            }

            for (int i = 0; i < ds.Count; i++)
            {
                j = (i - 1 + m_DomainPolygon.Count) % m_DomainPolygon.Count;
                k = (i + 1) % m_DomainPolygon.Count;

                if (ds[j] == 0)
                    ss.Add(0);
                else
                    ss.Add(ds[j] / (ds[j] + ds[k]));
            }

            for (int i = 0; i < ds.Count; i++)
                Distance.Add((ss[i], ds[i]));

            return Distance;
        }

        private List<(double, double)> ComputeRadialDistanceFunctionforKato(double u, double v, int n)
        {
            List<(double, double)> output = new List<(double, double)>();
            int i_before = 0, i_after = 0, i_afterafter = 0;
            int j;

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

                j = (i + 1) % (m_DomainPolygon.Count);
                //output.Add(((new Point3d(u, v, 0) - e_i).Length, (e_i - m_DomainPolygon[i]).Length));
                output.Add(((e_i - m_DomainPolygon[i]).Length / (m_DomainPolygon[j] - m_DomainPolygon[i]).Length, (new Point3d(u, v, 0) - e_i).Length));
            }
            return output;
        }

        private static void ComputeDomainPolygon()
        {
            /*double TotL = 0, L;
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
            }*/

            m_DomainPolygon = new List<Point3d>();
            for (int i = 0; i < m_Curves.Count; i++)
                m_DomainPolygon.Add(new Point3d(m_Curves[i].PointAtStart.X, m_Curves[i].PointAtStart.Y, 0)); //m_DomainPolygon.Add(new Point3d(0.1 * m_Curves[i].PointAtStart.X, 0.1 * m_Curves[i].PointAtStart.Y, 0));

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

        private void writerfunction(HarmonicMap map)
        {
            var filename = "C:\\Users\\SCAD\\Downloads\\file.csv";
            StreamWriter f = new StreamWriter(filename);
            for (double u = 0; u < 6.05; u = u + 0.2)
            {
                for (double v = 0; v < 6.05; v = v + 0.2)
                {
                    Point3d point2 = new Point3d(u, v, 0);
                    var success = harmonic_eval(map, point2, out double result);
                    if (success)
                    {
                        //f.Write("{0:N6},{0:N6},{0:N6}\n", u, v, result);
                        f.Write(u.ToString("N6") + "," + v.ToString("N6") + "," + result.ToString("N6") + '\n');
                        //f.Write(u.ToString() + ',' + v.ToString() + ',' + result.ToString() + '\n');
                    }
                }
            }
            f.Close();
        }
        private string writeroutput(HarmonicMap map)
        {
            //var filename = "C:\\Users\\SCAD\\Downloads\\file.csv";
            string f = "";
            for (double u = 0; u < 10.1; u = u + 0.5)
            {
                for (double v = 0; v < 20.1; v = v + 0.5)
                {
                    Point3d point2 = new Point3d(u, v, 0);
                    var success = harmonic_eval(map, point2, out double result);
                    if (success)
                    {
                        //f.Write("{0:N6},{0:N6},{0:N6}\n", u, v, result);
                        f += (u.ToString("N6") + "," + v.ToString("N6") + "," + result.ToString("N6") + '\n');
                        //f.Write(u.ToString() + ',' + v.ToString() + ',' + result.ToString() + '\n');
                    }
                }
            }
            return f;
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