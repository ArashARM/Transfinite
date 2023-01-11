using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;

namespace Grasshopper2
{
    public class Harmonic : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the Harmonic class.
        /// </summary>
        public Harmonic()
          : base("Harmonic", "Harmonic",
              "Harmonic",
              "Category", "Harmonic")
        {
        }

        private static List<Point3d> m_Pts;
        private static List<double> m_Values;

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddTextParameter("input", "input", "input",GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddPointParameter("Pts", "Pts", "Pts", GH_ParamAccess.list);
            pManager.AddNumberParameter("Values", "Values", "Values", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            string readdata = "";
            if (!DA.GetData(0,ref readdata)) return;
            if(readdata == "") ReadFile("");
            else ReadInput(readdata);

            double MinVal = 9E9, MaxVal = -9E9;
            for (int i = 0; i < m_Values.Count; i++)
            {
                if (m_Values[i] < MinVal)
                    MinVal = m_Values[i];
                if (m_Values[i] > MaxVal)
                    MaxVal = m_Values[i];
            }

            //MinVal = 0.001;
            //MaxVal = 1;

            for (int i = 0; i < m_Values.Count; i++)
                m_Values[i] = Math.Abs(MaxVal - (m_Values[i] - MinVal)) / MaxVal * 255;

            DA.SetDataList(0, m_Pts);
            DA.SetDataList(1, m_Values);
        }

        private static void ReadFile(string input)
        {
            m_Pts = new List<Point3d>();
            m_Values = new List<double>();
            string EAs;
            if (input == "")
            {
                EAs = System.IO.File.ReadAllText("C://Downloads//file.csv");
            }
            else EAs = input;
            string[] lines = EAs.Split(new string[] { "\r\n" }, StringSplitOptions.None);
            string[] coords;
            for (int i = 0; i < lines.Length; i++)
            {
                System.Diagnostics.Debug.WriteLine(i);

                coords = lines[i].Split(new string[] { "," }, StringSplitOptions.None);
                if (coords.Length == 1)
                    break;

                if (Convert.ToDouble(coords[0]) > 3.01 && Convert.ToDouble(coords[1]) > 3.01)
                    continue;

                m_Pts.Add(new Point3d(Convert.ToDouble(coords[0]), Convert.ToDouble(coords[1]), 0));
                m_Values.Add(Convert.ToDouble(coords[2]));               
            }
        }
        private static void ReadInput(string input)
        {
            m_Pts = new List<Point3d>();
            m_Values = new List<double>();
            string EAs = input;
            string[] lines = EAs.Split(new string[] { "\n" }, StringSplitOptions.None);
            string[] coords;
            for (int i = 0; i < lines.Length; i++)
            {
                System.Diagnostics.Debug.WriteLine(i);

                coords = lines[i].Split(new string[] { "," }, StringSplitOptions.None);
                if (coords.Length == 1)
                    break;

                if (Convert.ToDouble(coords[0]) > 5 && Convert.ToDouble(coords[1]) > 10)
                    continue;

                m_Pts.Add(new Point3d(Convert.ToDouble(coords[0]), Convert.ToDouble(coords[1]), 0));
                m_Values.Add(Convert.ToDouble(coords[2]));
            }
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
            get { return new Guid("36f400de-6fdf-4688-b6f3-0150b11d932b"); }
        }
    }
}