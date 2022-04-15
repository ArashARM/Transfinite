using Grasshopper;
using Grasshopper.Kernel.Data;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.Serialization.Formatters.Binary;
using System.Text;
using System.Threading.Tasks;
using Accord.Math;

namespace Extensions
{
    public static class SurfaceExtensions
    {
        public static void Reparameterize(this NurbsSurface surface)
        {
            var dm = new Interval(0, 1);
            surface.SetDomain(0, dm);
            surface.SetDomain(1, dm);
        }
        public static void Reparameterize(this NurbsCurve curve)
        {
            var dm = new Interval(0, 1);
            curve.Domain = dm;
        }

        public static void Reparameterize(this Curve curve)
        {
            var dm = new Interval(0, 1);
            curve.Domain = dm;
        }

        //public static Curve Reparameterize(this Curve curve)
        //{
        //    var dm = new Interval(0, 1);
        //    curve.Domain = dm;
        //    return curve;
        //}

        public static void Reparameterize(this IEnumerable<NurbsCurve> curve)
        {
            var dm = new Interval(0, 1);
            foreach (var item in curve)
                item.Domain = dm;
        }
        public static void Reparameterize(this IEnumerable<Curve> curve)
        {
            var dm = new Interval(0, 1);
            foreach (var item in curve)
                item.Domain = dm;
        }

        public static NurbsSurface TransposeToNurbs(this NurbsSurface surface)
        {
            Surface newsurf = surface;
            surface = newsurf.Transpose().ToNurbsSurface();
            return surface;
        }

        public static NurbsSurface ReverseToNurbs(this NurbsSurface surface, int direction)
        {
            Surface newsurf = surface;

            if (direction == 0)
                surface = newsurf.Reverse(0).ToNurbsSurface();
            else
                surface = newsurf.Reverse(1).ToNurbsSurface();

            return surface;

        }

        public static void EqualizeKnots(ref NurbsSurface A, ref NurbsSurface B)
        {
            A.Reparameterize();
            B.Reparameterize();

            var distknotsA_U = A.KnotsU.DistinctKnots(B.KnotsU);
            var distknotsB_U = B.KnotsU.DistinctKnots(A.KnotsU);

            var distknotsA_V = A.KnotsV.DistinctKnots(B.KnotsV);
            var distknotsB_V = B.KnotsV.DistinctKnots(A.KnotsV);

            foreach (var item in distknotsA_U)
                B.KnotsU.InsertKnot(item);
            foreach (var item in distknotsB_U)
                A.KnotsU.InsertKnot(item);
            foreach (var item in distknotsA_V)
                B.KnotsV.InsertKnot(item);
            foreach (var item in distknotsB_V)
                A.KnotsV.InsertKnot(item);
        }

        public static void EqualizeKnots(ref Surface ASurf, ref Surface BSurf)
        {

            var A = ASurf.ToNurbsSurface();
            A.Reparameterize();
            var B = BSurf.ToNurbsSurface();
            B.Reparameterize();
            var distknotsA_U = A.KnotsU.DistinctKnots(B.KnotsU);
            var distknotsB_U = B.KnotsU.DistinctKnots(A.KnotsU);

            var distknotsA_V = A.KnotsV.DistinctKnots(B.KnotsV);
            var distknotsB_V = B.KnotsV.DistinctKnots(A.KnotsV);

            foreach (var item in distknotsA_U)
                B.KnotsU.InsertKnot(item);
            foreach (var item in distknotsB_U)
                A.KnotsU.InsertKnot(item);
            foreach (var item in distknotsA_V)
                B.KnotsV.InsertKnot(item);
            foreach (var item in distknotsB_V)
                A.KnotsV.InsertKnot(item);

            ASurf = A;
            BSurf = B;
        }

    }

    public static class LineExtensions
    {
        public static double ParalelPoint(this Rhino.Geometry.Line line)
        {
            if (line.FromX == line.ToX)
                return line.ToX;
            else if (line.FromY == line.ToY)
                return line.ToY;
            else
                return 0;
        }

        public static Rhino.Geometry.Line MovetoSide(this Rhino.Geometry.Line line, double t, bool left)
        {
            Transform move = new Transform();
            if (left)
                t = -t;
            Rhino.Geometry.Line movedline = line;
            if (line.FromX == line.ToX)
                move = Transform.Translation(t, 0, 0);
            else if (line.FromY == line.ToY)
                move = Transform.Translation(0, t, 0);

            movedline.Transform(move);
            return movedline;
        }

        public static Point3d MidPoint(this Rhino.Geometry.Line line)
        {
            return line.PointAt(0.5);
        }

        public static Point2d MidPointUV(this Rhino.Geometry.Line line)
        {
            return line.PointAtUV(0.5);
        }

        public static double Slope(this Rhino.Geometry.Line line)
        {
            double m;
            //if (line.ToX - line.FromX != 0)
            m = (line.ToY - line.FromY) / (line.ToX - line.FromX);

            m = Math.Atan(m);
            //else
            //{
            //    if (line.ToY - line.FromY < 0)
            //        m = double.MinValue;
            //    else
            //        m = double.MaxValue;
            //}
            return m;
        }

        public static Point2d PointAtUV(this Rhino.Geometry.Line line, double t)
        {
            Point2d pt;
            if (line.FromX == line.ToX)
                pt = new Point2d(line.From.X, line.PointAt(t).Y);
            else if (line.FromY == line.ToY)
                pt = new Point2d(line.PointAt(t).X, line.FromY);
            else
                pt = new Point2d(line.PointAt(t).X, line.PointAt(t).Y);

            return pt;
        }

        public static Rhino.Geometry.Line Extend(this Rhino.Geometry.Line line)
        {
            Line ln;
            ln = line;
            if (ln.FromX == ln.ToX)
            {
                if (ln.FromY > ln.ToY)
                {
                    ln.FromY = 1;
                    ln.ToY = 0;
                }
                else
                {
                    ln.FromY = 0;
                    ln.ToY = 1;
                }
            }
            else if (ln.FromY == ln.ToY)
            {
                if (ln.FromX > ln.ToX)
                {
                    ln.FromX = 1;
                    ln.ToX = 0;
                }
                else
                {
                    ln.FromX = 0;
                    ln.ToX = 1;
                }
            }


            return ln;
        }

    }

    public static class EnumExtensions
    {
        public static T SecondLast<T>(this IEnumerable<T> list)
        {
            var lst = list.ToList();
            var secondLast = lst[lst.Count() - 2];
            return secondLast;
        }

        //public static T SecondLast<T>(this List<T> list)
        //{
        //    var secondLast = list[list.Count() - 2];
        //    return secondLast;
        //}

        public static void AddToFront<T>(this List<T> list, T item)
        {
            // omits validation, etc.
            list.Insert(0, item);
        }

        public static void AddCurveEnds(this List<int> order, bool First)
        {
            if (order.Count < 1)
                return;
            if (First)
            {
                //first curve
                if (order.First<int>() % 2 == 0 && !order.Contains(order.First<int>() + 1))
                    order.AddToFront(order.First<int>() + 1);
                else if (order.First<int>() % 2 == 1 && !order.Contains(order.First<int>() - 1))
                    order.AddToFront(order.First<int>() - 1);
            }
            else
            {
                order.RemoveAt(0);
            }

            //LastCurve
            if (order.Last<int>() % 2 == 0 && !order.Contains(order.Last<int>() + 1))
                order.Add(order.Last<int>() + 1);
            else if (order.Last<int>() % 2 == 1 && !order.Contains(order.Last<int>() - 1))
                order.Add(order.Last<int>() - 1);

        }

    }

    public static class FileExtensions
    {
        /// <summary>
        /// Line Count for loops
        /// </summary>
        /// <param name="path"></param>
        /// <returns> it returns count -1 </returns>
        public static int FileLineCount(this string path)
        {
            if (File.Exists(path))
                return File.ReadLines(path).Count() - 1;
            else
                return 0;
        }
        public static string FileLastLine(this string path)
        {
            if (File.Exists(path))
                return File.ReadLines(path).Last();
            else
                return null;
        }

        /// <summary>
        /// Gets last line of path file
        /// </summary>
        /// <param name="path"></param>
        /// <param name="spr"></param>
        /// <returns> splited last line</returns>
        public static string[] FileLastLineSep(this string path, string spr = ";")
        {
            if (File.Exists(path))
                return File.ReadLines(path).Last().Split(spr.ToCharArray());
            else
                return null;
        }

        public static bool ChangeLine(this string filepath, string newline, int line_to_edit)
        {
            if (line_to_edit < 0)
            {

                string[] allLine = File.ReadAllLines(filepath);
                allLine[line_to_edit - 1] = newline;
                File.WriteAllLines(filepath, allLine);
                return true;
            }
            else
                return false;
        }

        public static void ChangeLastLine(this string filepath, string newline)
        {
            var line_to_edit = File.ReadLines(filepath).Count();
            string[] allLine = File.ReadAllLines(filepath);
            allLine[line_to_edit - 1] = newline;
            File.WriteAllLines(filepath, allLine);
        }

        //public static void RewriteData(this string filepath, IEnumerable<string> newlines)
        //{

        //    File.WriteAllLines(filepath, newlines);
        //}
    }

    public static class KnotExtensions
    {
        public static List<double> ToList(this Rhino.Geometry.Collections.NurbsSurfaceKnotList knots)
        {
            var knotList = new List<double>();
            foreach (var item in knots)
            {
                knotList.Add(item);
            }
            return knotList;
        }

        public static List<double> DistinctKnots(this Rhino.Geometry.Collections.NurbsSurfaceKnotList knots, Rhino.Geometry.Collections.NurbsSurfaceKnotList knots2)
        {

            var list1 = knots.ToList();
            var list2 = knots2.ToList();
            return list1.Except(list2).ToList();
        }


    }

    public static class PointExtensions
    {
        public static Point2d ToPoint2d(this Point3d point)
        {
            Point2d newpoint = new Point2d(point.X, point.Y);
            return newpoint;
        }

        public static Point3d ToPoint3d(this Point2d point)
        {
            Point3d newpoint = new Point3d(point.X, point.Y, 0);
            return newpoint;
        }

        public static List<Point2d> ToPoint2d(this IEnumerable<Point3d> point)
        {
            List<Point2d> newpoint = new List<Point2d>();
            foreach (var item in point)
            {
                newpoint.Add(new Point2d(item.X, item.Y));
            }
            return newpoint;
        }

        public static List<Point3d> ToPoint3d(this IEnumerable<Point2d> point)
        {
            List<Point3d> newpoint = new List<Point3d>();
            foreach (var item in point)
            {
                newpoint.Add(new Point3d(item.X, item.Y, 0));
            }
            return newpoint;
        }

        public static Point3d MoveX(this Point3d point, double translation_x)
        {
            var pt = new Point3d(point.X + translation_x, point.Y, point.Z);
            return pt;
        }

        public static Point3d MoveY(this Point3d point, double translation_y)
        {
            var pt = new Point3d(point.X, point.Y + translation_y, point.Z);
            return pt;
        }

        public static Point3d MoveZ(this Point3d point, double translation_z)
        {
            var pt = new Point3d(point.X, point.Y, point.Z + translation_z);
            return pt;
        }

    }

    public static class ListExtensions
    {
        public static void Mirror<T>(this T[,] array, int i, int j)
        {
            array[j, i] = array[i, j];
        }
        public static void Mirror<T>(this T[,] array)
        {
            for (int i = 0; i < array.Rows(); i++)
            {
                for (int j = i; j < array.Columns(); j++)
                {
                    array[j, i] = array[i, j];
                }
            }
        }
        public static T DeepClone<T>(this T obj)
        {
            using (var ms = new MemoryStream())
            {
                var formatter = new BinaryFormatter();
                formatter.Serialize(ms, obj);
                ms.Position = 0;
                return (T)formatter.Deserialize(ms);
            }
        }


        public static DataTree<T> ToTree<T>(this List<List<T>> list)
        {
            var tree = new DataTree<T>();

            for (int i = 0; i < list.Count; i++)
            {
                var pth = new GH_Path(i);
                foreach (var result in list[i])
                {
                    tree.Add(result, pth);
                }
            }
            return tree;
        }

        public static void RemoveAt(this List<int> lst, IEnumerable<int> extractIdxs)
        {
            //var ordered = extractIdxs.OrderByDescending(i => i);
            foreach (var item in extractIdxs)
            {
                if (item % 2 == 0 && !extractIdxs.Contains(item + 1))
                    continue;
                else if (item % 2 == 1 && !extractIdxs.Contains(item - 1))
                    continue;
                else
                    lst.Remove(item);
            }
        }
    }

    public static class DoubleExtensions
    {
        public static string ToS(this double num, int dec = 3)
        {
            string deci = "0";
            for (int i = 0; i < dec; i++)
            {
                if (i == 0)
                    deci += ".";
                deci += "0";
            }

            return num.ToString(deci);
        }
    }
    public static class IntExtensions
    {
        public static string ToS(this int num, int dec = 3)
        {
            string deci = "0";
            for (int i = 0; i < dec; i++)
            {
                if (i == 0)
                    deci += ".";
                deci += "0";
            }

            return num.ToString(deci);
        }
    }
}


