using System;
using System.Collections.Generic;
using System.IO;
using Grasshopper.Kernel;
using Rhino.Geometry;
using static GrasshopperProjects.Harmonic;

namespace GrasshopperProjects
{
    public class HarmonicComponent : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the Harmonic class.
        /// </summary>
        public HarmonicComponent()
          : base("Harmonic", "Nickname",
              "Description",
              "GHP", "Transfinite")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            /*double input[] = { 0, 0, 0,
		    6, 0, 0,
		    6, 6, 0,
		    3.5, 6, 0,
		    3, 3, 1,
		    2.5, 6, 0,
		    0, 6, 0 };*/
            Point3d[] input = {   new Point3d(0, 0, 1),
                                  new Point3d(6, 0, 0),
                                  new Point3d(6, 3, 0),
                                  new Point3d(3, 3, 0),
                                  new Point3d(3, 6, 1),
                                  new Point3d(0, 6, 0)};

            double[] min = { 0, 0 }, max = { 6, 6 };
            bool success;
            double result;
            Point3d point = new Point3d(2, 2, 0);
            int levels = 9;
            HarmonicMap map;
            double u, v;

            map = harmonic_create(min, max, levels);

            /*for (int i = 0; i < 7; ++i) {
                int j = (i + 1) % 7;
                harmonic_add_point(map, &input[3 * i]);
            }*/

            for (int i = 0; i < 6; ++i)
            {
                int j = (i + 1) % 6;
                harmonic_add_line(map, input[i], input[j]);
            }
            /*for (int i = 0; i < 6; ++i) {
                int j = (i + 1) % 6;
                harmonic_add_line(map, &input[3 * i], &input[3 * j]);
            }*/
            harmonic_solve(map, 1.0e-5, false);

            /* Evaluation test */
            success = harmonic_eval(map, point, out result);
            //if (success)
            //    //printf("%lf\n", result);
            //else
            //    printf("N/A\n");
            var filename = "C://Users//alper//Desktop//Hrmnc3//Test//Test//file.csv";
            StreamWriter f = new StreamWriter(filename);
            for (u = 0; u < 6.05; u = u + 0.2)
            {
                for (v = 0; v < 6.05; v = v + 0.2)
                {
                    Point3d point2 = new Point3d(u, v, 0);
                    success = harmonic_eval(map, point2, out result);
                    if (success)
                    {
                        //f.Write("{0:N6},{0:N6},{0:N6}\n", u, v, result);
                        f.Write(u.ToString("N6") + "," + v.ToString("N6") + "," + result.ToString("N6") + '\n');
                        //f.Write(u.ToString() + ',' + v.ToString() + ',' + result.ToString() + '\n');
                    }
                }
            }
            f.Close();
            /* PPM output test */
            harmonic_write_ppm(map, "C://Users//alper//Desktop//Hrmnc3//Test//Test//test.ppm");

            harmonic_free(map);
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
            get { return new Guid("A6E2851D-0409-46AC-8809-B356981FD13B"); }
        }
    }
}