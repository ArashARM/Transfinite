using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Rhino.Geometry;

namespace GrasshopperProjects
{
    public class Harmonic
    {
        static int MIN_LEVEL = 3;
        static double MAX(double a, double b) => (a > b) ? a : b;

        public struct GridValue
        {
            public bool boundary { get; set; }
            public double value { get; set; }
        };

        public struct HarmonicMap
        {
            public int levels { get; set; }
            public int size { get; set; }
            public GridValue[] grid { get; set; }
            public double[] offset { get; set; }
            public double scaling { get; set; }

        };


        public static HarmonicMap harmonic_create(double[] min, double[] max, int levels)
        {
            /* Create the grid */
            int n = (int)Math.Pow(2, levels);
            HarmonicMap map = new HarmonicMap();
            map.offset = new double[2];
            map.levels = levels;
            map.size = n;
            map.grid = new GridValue[n * n];

            /* Add a margin of 2.5% on all sides */
            double length = MAX(max[0] - min[0], max[1] - min[1]);
            map.offset[0] = min[0] - length * 0.025;
            map.offset[1] = min[1] - length * 0.025;
            map.scaling = (double)n / length / 1.05;

            /* Initialize cells */
            for (int i = 0; i < n * n; ++i)
            {
                map.grid[i].boundary = false;
                map.grid[i].value = 0.0;
            }
            return map;
        }

        public static void solveHarmonic(GridValue[] grid, int n, double epsilon)
        {
            double change;
            do
            {
                change = 0.0;
                int count = 0, index = n + 1;
                for (int j = 1, n_1 = n - 1; j < n_1; ++j)
                {
                    for (int i = 1; i < n_1; ++i, ++index)
                        if (!grid[index].boundary)
                        {
                            double value = 0.0;
                            value += grid[index - n].value;
                            value += grid[index - 1].value;
                            value += grid[index + n].value;
                            value += grid[index + 1].value;
                            value /= 4.0;
                            change += Math.Abs(grid[index].value - value);
                            grid[index].value = value;
                            ++count;
                        }
                    index += 2;
                }
                /* Boundary cases: not handled [the result is the same] */
                change /= (double)count;
            } while (change > epsilon);
        }

        public static void solveBiharmonic(GridValue[] grid, int n, double epsilon)
        {
            double change;
            do
            {
                change = 0.0;
                int count = 0, n_2 = n - 2, n2 = n * 2, index = n2 + 2;
                for (int j = 2; j < n_2; ++j)
                {
                    for (int i = 2; i < n_2; ++i, ++index)
                        if (!grid[index].boundary)
                        {
                            double value = 0.0;
                            value += grid[index - n].value;
                            value += grid[index - 1].value;
                            value += grid[index + n].value;
                            value += grid[index + 1].value;
                            value *= 4.0;
                            value -= grid[index - n - 1].value;
                            value -= grid[index - n + 1].value;
                            value -= grid[index + n - 1].value;
                            value -= grid[index + n + 1].value;
                            value *= 2.0;
                            value -= grid[index - n2].value;
                            value -= grid[index - 2].value;
                            value -= grid[index + n2].value;
                            value -= grid[index + 2].value;
                            value /= 20.0;
                            change += Math.Abs(grid[index].value - value);
                            grid[index].value = value;
                            ++count;
                        }
                    index += 4;
                }
                /* Boundary cases */
                for (int j = 0; j < n; ++j)
                    for (int i = 0; i < n; ++i)
                    {
                        if (j >= 2 && j < n_2 && i >= 2 && i < n_2)
                            continue;
                        index = j * n + i;
                        if (grid[index].boundary)
                            continue;
                        double weight = 0.0;
                        double value = 0.0;
                        if (i > 0)
                        {
                            int k = index - 1;
                            weight += 2.0;
                            value += grid[k].value * 2.0;
                            int neighbors = 4;
                            if (i == 1)
                                --neighbors;
                            if (j == 0)
                                --neighbors;
                            if (j == n - 1)
                                --neighbors;
                            double w = 1.0 / (double)neighbors;
                            if (i > 1)
                            {
                                weight -= w;
                                value -= grid[k - 1].value * w;
                            }
                            if (j > 0)
                            {
                                weight -= w;
                                value -= grid[k - n].value * w;
                            }
                            if (j < n - 1)
                            {
                                weight -= w;
                                value -= grid[k + n].value * w;
                            }
                        }
                        if (i < n - 1)
                        {
                            int k = index + 1;
                            weight += 2.0;
                            value += grid[k].value * 2.0;
                            int neighbors = 4;
                            if (i == n - 2)
                                --neighbors;
                            if (j == 0)
                                --neighbors;
                            if (j == n - 1)
                                --neighbors;
                            double w = 1.0 / (double)neighbors;
                            if (i < n - 2)
                            {
                                weight -= w;
                                value -= grid[k + 1].value * w;
                            }
                            if (j > 0)
                            {
                                weight -= w;
                                value -= grid[k - n].value * w;
                            }
                            if (j < n - 1)
                            {
                                weight -= w;
                                value -= grid[k + n].value * w;
                            }
                        }
                        if (j > 0)
                        {
                            int k = index - n;
                            weight += 2.0;
                            value += grid[k].value * 2.0;
                            int neighbors = 4;
                            if (i == 0)
                                --neighbors;
                            if (i == n - 1)
                                --neighbors;
                            if (j == 1)
                                --neighbors;
                            double w = 1.0 / (double)neighbors;
                            if (i > 0)
                            {
                                weight -= w;
                                value -= grid[k - 1].value * w;
                            }
                            if (i < n - 1)
                            {
                                weight -= w;
                                value -= grid[k + 1].value * w;
                            }
                            if (j > 1)
                            {
                                weight -= w;
                                value -= grid[k - n].value * w;
                            }
                        }
                        if (j < n - 1)
                        {
                            int k = index + n;
                            weight += 2.0;
                            value += grid[k].value * 2.0;
                            int neighbors = 4;
                            if (i == 0)
                                --neighbors;
                            if (i == n - 1)
                                --neighbors;
                            if (j == n - 2)
                                --neighbors;
                            double w = 1.0 / (double)neighbors;
                            if (i > 0)
                            {
                                weight -= w;
                                value -= grid[k - 1].value * w;
                            }
                            if (i < n - 1)
                            {
                                weight -= w;
                                value -= grid[k + 1].value * w;
                            }
                            if (j < n - 2)
                            {
                                weight -= w;
                                value -= grid[k + n].value * w;
                            }
                        }
                        value /= weight;
                        change += Math.Abs(grid[index].value - value);
                        grid[index].value = value;
                        ++count;
                    }
                change /= (double)count;
            } while (change > epsilon);
        }

        public static void solve(GridValue[] grid, int level, double epsilon, bool biharmonic)
        {
            int n = (int)Math.Pow(2, level);
            if (level > MIN_LEVEL)
            {
                /* Generate a coarser grid and solve that first to get good starting values */
                int level1 = level - 1, n1 = (int)Math.Pow(2, level1);
                GridValue[] grid1 = new GridValue[n1 * n1];
                for (int i = 0; i < n1; ++i)
                    for (int j = 0; j < n1; ++j)
                    {
                        grid1[j * n1 + i].value = 0.0;
                        int count = 0;
                        if (grid[2 * j * n + 2 * i].boundary)
                        {
                            ++count;
                            grid1[j * n1 + i].value += grid[2 * j * n + 2 * i].value;
                        }
                        if (grid[2 * j * n + 2 * i + 1].boundary)
                        {
                            ++count;
                            grid1[j * n1 + i].value += grid[2 * j * n + 2 * i + 1].value;
                        }
                        if (grid[(2 * j + 1) * n + 2 * i].boundary)
                        {
                            ++count;
                            grid1[j * n1 + i].value += grid[(2 * j + 1) * n + 2 * i].value;
                        }
                        if (grid[(2 * j + 1) * n + 2 * i + 1].boundary)
                        {
                            ++count;
                            grid1[j * n1 + i].value += grid[(2 * j + 1) * n + 2 * i + 1].value;
                        }
                        if (count > 0)
                        {
                            grid1[j * n1 + i].boundary = true;
                            grid1[j * n1 + i].value /= (double)count;
                        }
                        else
                            grid1[j * n1 + i].boundary = false;
                    }
                solve(grid1, level1, epsilon, biharmonic);
                for (int i = 0; i < n; ++i)
                    for (int j = 0; j < n; ++j)
                        if (!grid[j * n + i].boundary)
                            grid[j * n + i].value = grid1[(j / 2) * n1 + i / 2].value;
                grid1 = null;
            }

            /* Solve by iteration */
            if (biharmonic)
                solveBiharmonic(grid, n, epsilon);
            else
                solveHarmonic(grid, n, epsilon);
        }

        public static void harmonic_add_point(HarmonicMap map, Point3d point)
        {
            int n = map.size;
            int x = (int)Math.Round((point[0] - map.offset[0]) * map.scaling);
            int y = (int)Math.Round((point[1] - map.offset[1]) * map.scaling);
            map.grid[y * n + x].boundary = true;
            map.grid[y * n + x].value = point[2];
        }

        public static void harmonic_add_line(HarmonicMap map, Point3d from, Point3d to)
        {
            int n = map.size;
            int x0 = (int)Math.Round((from[0] - map.offset[0]) * map.scaling);
            int y0 = (int)Math.Round((from[1] - map.offset[1]) * map.scaling);
            double v0 = from[2];
            int x1 = (int)Math.Round((to[0] - map.offset[0]) * map.scaling);
            int y1 = (int)Math.Round((to[1] - map.offset[1]) * map.scaling);
            double v1 = to[2];
            int dx = Math.Abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
            int dy = Math.Abs(y1 - y0), sy = y0 < y1 ? 1 : -1;
            int err = (dx > dy ? dx : -dy) / 2, e2;
            if (err == 0)
            {
                map.grid[y0 * n + x0].boundary = true;
                map.grid[y0 * n + x0].value = v0;
                map.grid[y1 * n + x1].boundary = true;
                map.grid[y1 * n + x1].value = v1;
                return;
            }
            while (true)
            {
                double ratio;             /* linear interpolation along the sides */
                if (err > 0)
                    ratio = (double)Math.Abs(x1 - x0) / (double)dx;
                else
                    ratio = (double)Math.Abs(y1 - y0) / (double)dy;
                map.grid[y0 * n + x0].boundary = true;
                map.grid[y0 * n + x0].value = v0 * ratio + v1 * (1.0 - ratio);
                if (x0 == x1 && y0 == y1) break;
                e2 = err;
                if (e2 > -dx) { err -= dy; x0 += sx; }
                if (e2 < dy) { err += dx; y0 += sy; }
            }
        }

        public static void harmonic_add_curve(HarmonicMap map, Point3d[] points, int n)
        {
            double tmp;
            int resolution = map.size;
            double[] coeff = new double[n];
            Point3d from, to;
            from = new Point3d(points[0].X, points[0].Y, points[0].Z);
            for (int i = 1; i <= resolution; ++i)
            {
                double u = (double)i / resolution;
                /* Compute Bernstein polynomials */
                coeff[0] = 1.0;
                for (int j = 1; j < n; ++j)
                {
                    double saved = 0.0;
                    for (int k = 0; k < j; ++k)
                    {
                        tmp = coeff[k];
                        coeff[k] = saved + tmp * (1.0 - u);
                        saved = tmp * u;
                    }
                    coeff[j] = saved;
                }
                /* Evaluate the curve */
                to = new Point3d(0.0, 0.0, 0.0);
                for (int j = 0; j < n; ++j)
                {
                    to += points[j] * coeff[j];
                }
                /* Draw a segment */
                harmonic_add_line(map, from, to);
                /* Swap from & to */
                (from, _) = (to, from);
            }
            coeff = null;
        }

        public static void harmonic_solve(HarmonicMap map, double epsilon, bool biharmonic)
        {
            solve(map.grid, map.levels, epsilon, biharmonic);
        }

        public static bool inside_map(HarmonicMap map, int i, int j)
        {
            return i >= 0 && j >= 0 && i < map.size && j < map.size;
        }

        public static bool harmonic_eval(HarmonicMap map, Point3d point, out double value)
        {
            int n = map.size;
            double x = (point[0] - map.offset[0]) * map.scaling;
            double y = (point[1] - map.offset[1]) * map.scaling;
            int i = (int)Math.Round(x), j = (int)Math.Round(y);

            if (!(inside_map(map, i, j) &&
                  inside_map(map, i, j + 1) &&
                  inside_map(map, i + 1, j) &&
                  inside_map(map, i + 1, j + 1)))
            {
                value = 0;
                return false;               /* The point is outside the region */
            }

            value = map.grid[j * n + i].value * (1.0 - y + j) * (1.0 - x + i);
            value += map.grid[(j + 1) * n + i].value * (y - j) * (1.0 - x + i);
            value += map.grid[j * n + i + 1].value * (1.0 - y + j) * (x - i);
            value += map.grid[(j + 1) * n + i + 1].value * (y - j) * (x - i);

            return true;
        }

        public static void harmonic_write_ppm(HarmonicMap map, string filename)
        {
            int n = map.size;
            StreamWriter f = new StreamWriter(filename);

            f.Write("P3\n%zu %zu\n255\n", n, n);
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < n; ++j)
                    if (map.grid[j * n + i].boundary)
                        f.Write("255 0 0 ");
                    else
                        f.Write("0 0 " + ((int)Math.Round(map.grid[j * n + i].value * 255.0)).ToString("N6"));
                //f.Write("0 0 %d ", (int)Math.Round(map.grid[j * n + i].value * 255.0));
                f.Write("\n");
            }
            f.Close();
        }

        public static void harmonic_free(HarmonicMap map)
        {
            map.grid = null;
            map = new HarmonicMap();
        }

    }
}
