using ComplexOperations;
using System;
using System.IO;

namespace Tunneling_Effect_2._0
{
    class Program
    {
        static Complex AC, Bconst; // parameters of Cranc Nichols methods
        static Complex[] B;
        static Complex[] DPotential;
        static void Main(string[] args)
        {
            string path = @"C:\Users\Алеша\Documents\Tunneling_Effect_2.0";
            int N = 40000;
            double waveNumber = 7.07, waveWidth = 3, potentialWidth = 3, PotentialMaxValue = 25;
            double waveStartPoint = 40 + 4 * waveWidth;
            double potentialCenter = 40 + 8 * waveWidth + 4 * potentialWidth + 10;
            double dx, dt, waveLength;
            double[] param = new double[5] { 7, 32, 12, 32, 14 };

            waveLength = 2 * Math.PI / waveNumber;
            dx = waveLength / 100;
            dt = dx / (2 * waveNumber);

            Console.WriteLine($"dx={dx}, dt={dt}");

            Complex i = new Complex(0, 1);
            Complex[] waveFunction = new Complex[N + 1];
            B = new Complex[N + 1];
            DPotential = new Complex[N + 1];

            AC = i * dt / 4 / dx / dx;
            Bconst = 1 + 2 * AC;


            using (StreamWriter sw = new StreamWriter($"{path}/Potential.dat", false, System.Text.Encoding.Default))
            {
                for (int j = 0; j < N + 1; j++)
                {
                    double xCoordinate = j * dx;
                    double potential = PotentialMaxValue * Math.Exp(
                                    -(Math.Pow((xCoordinate - potentialCenter) / potentialWidth, 2)
                                    ));


                    waveFunction[j] = Complex.Exp(i * waveNumber * xCoordinate) *
                                      Math.Exp(
                                          -(Math.Pow((xCoordinate - waveStartPoint) / waveWidth, 2)
                                          ));

                    DPotential[j] = potential * i * dt / 2;
                    B[j] = Bconst + DPotential[j];
                    DPotential[j] = 1 - DPotential[j];

                    sw.WriteLine($"{xCoordinate.ToString().Replace(',', '.')}"
                        + " " + $"{potential.ToString().Replace(',', '.')}");
                }
            }

                double Summa = Sum(N, waveFunction, param, dx);

                Console.WriteLine($"Summa = { Summa } ");
                Console.WriteLine(Energy(N, dx, waveFunction, param) / Summa);

                int Ntime = 0;
                int Npri = 200;
                int Mcicle = 0;

                while (true)
                {
                    Step(N, waveFunction);
                    Ntime++;
                    if (Ntime < Npri) continue;

                    Mcicle++;
                    Ntime = 0;

                    Console.WriteLine($"Mcicle = {Mcicle}");
                    Console.WriteLine($"Integral = {Sum(N, waveFunction, param, dx)}");
                    Console.WriteLine(Energy(N, dx, waveFunction, param) / Summa);
                    Console.WriteLine();


                using (StreamWriter sw = new StreamWriter($"{path}/Distrib.dat", false, System.Text.Encoding.Default))
                {
                    for (int j = 0; j < N + 1; j++)
                    {
                        sw.WriteLine($"{(j*dx).ToString().Replace(',', '.')}" + " " +
                        $"{Complex.GetReal(waveFunction[j]).ToString().Replace(',','.')}" + " " +
                        $"{Complex.GetImag(waveFunction[j]).ToString().Replace(',', '.')}" + " " +
                        $"{waveFunction[j].Norm().ToString().Replace(',', '.')}");
                    }
                }
                Console.ReadKey();
                }
            }

            private static double Energy(int N, double dx, Complex[] waveFunction, double[] param)
            {
                Complex Energy = new Complex(0, 0);
                int count = -1;
                Complex temp, temp1;

                for (int j = 1; j < N; j++)
                {
                    temp1 = waveFunction[j + 1] - 2 * waveFunction[j] + waveFunction[j - 1];
                    temp = waveFunction[j].Conj() * temp1;

                    count++;
                    if (count == 5) count = 1;

                    Complex add = temp * param[count];

                    Energy = Energy - add;
                }

                Energy = Energy / 45 / dx;

                return Complex.GetReal(Energy);
            }
            public static void Step(int N, Complex[] func)
            {
                Complex[] p = new Complex[N];
                Complex[] q = new Complex[N];
                Complex D, denominator;

                p[0] = new Complex(0, 0);
                q[0] = new Complex(0, 0);

                for (int count = 1; count < N; count++)
                {
                    D = DPotential[count] * func[count]
                        + AC * (func[count + 1] - 2 * func[count] + func[count - 1]);

                    denominator = B[count] - AC * p[count - 1];

                    p[count] = AC / denominator;
                    q[count] = (D + AC * q[count - 1]) / denominator;
                }

                func[N] = new Complex(0, 0);

                for (int count = N - 1; count >= 0; count--)
                {
                    func[count] = func[count + 1] * p[count] + q[count];
                }
            }
            public static double Sum(int N, Complex[] func, double[] param, double dx)
            {
                double Summa = 0;
                int count = -1;

                for (int j = 0; j < N + 1; j++)
                {
                    count++;
                    if (count == 5) count = 1;
                    Summa += func[j].Norm() * param[count];
                }

                Summa -= 7 * func[N].Norm();

                Summa = Summa * dx / 22.5;

                return Summa;
            }
        }
    }
