using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Numerics;

namespace FrameworkLibrary
{

        public interface Field<T>
        {
            int Dim { get; set; }
            double[] Size { get; set; }
            int[] N { get; set; }
            T Value { get; set; }
        }
        class Wavefunction : Field<Complex[,,]>
        {
            ////////////////////////
            //Property Declarations
            ///////////////////////


            private int dim;
            public int Dim
            {
                get
                {
                    return dim;
                }
                set
                {
                    dim = value;
                }

            }
            private double[] size;
            public double[] Size
            {
                get
                {
                    return size;
                }
                set
                {
                    size = value;
                }

            }
            private int[] n;
            public int[] N
            {
                get
                {
                    return n;
                }
                set
                {
                    n = value;
                }

            }
            private double[][] value;
            public double[][] Value
            {
                get
                {
                    return this.value;
                }
                set
                {
                    this.value = value;
                }

            }

            ////////////////////
            /// Private Methods
            ////////////////////

            private Wavefunction(int dimension, double[] window, int[] samplesize)
            {

                // Initialize field as zero everywhere

                int i;
                Dim = dimension;
                N = new int[dim];
                Size = new double[dim];
                for (i = 0; i < dim; i++)
                {
                    N[i] = samplesize[i];
                }
                for (i = 0; i < dim; i++)
                {
                    Value[i] = new double[N[i]];
                }
                Size = window;
                foreach (double[] f in Value)
                {
                    for (i = 0; i < dim; i++) { f[i] = 0; }
                }
            }



        }
        class Potential : Field<double[,,]>
        {
            ////////////////////////
            //Property Declarations
            ///////////////////////


            private int dim;
            public int Dim
            {
                get
                {
                    return dim;
                }
                set
                {
                    dim = value;
                }

            }
            private double[] size;
            public double[] Size
            {
                get
                {
                    return size;
                }
                set
                {
                    size = value;
                }

            }
            private int[] n;
            public int[] N
            {
                get
                {
                    return n;
                }
                set
                {
                    n = value;
                }

            }
            private double[,,] value;
            public double[,,] Value
            {
                get
                {
                    return this.value;
                }
                set
                {
                    this.value = value;
                }

            }

            ////////////////////
            /// Private Methods
            ////////////////////

            private Potential(int dimension, double[] window, int[] samplesize)
            {

                // Initialize field as zero everywhere

                int i;
                int j;
                int m;
                Dim = dimension;
                N = new int[dim];
                Size = new double[dim];
                N = samplesize;
                Size = window;
                Value = new double[N[0],N[1],N[2]];
                for (i = 0; i < N[0]; i++ )
                {
                    for (j = 0; j < N[0]; j++)
                    {
                        for (m = 0; m < N[0]; m++)
                        {
                            Value[i, j, m] = 0;
                        }
                    }
                }
                    

            }



        }
        class Electric : Field<double[][,,]>
        {
            ////////////////////////
            //Property Declarations
            ///////////////////////


            private int dim;
            public int Dim
            {
                get
                {
                    return dim;
                }
                set
                {
                    dim = value;
                }

            }
            private double[] size;
            public double[] Size
            {
                get
                {
                    return size;
                }
                set
                {
                    size = value;
                }

            }
            private int[] n;
            public int[] N
            {
                get
                {
                    return n;
                }
                set
                {
                    n = value;
                }

            }
            private double[][, ,] value;
            public double[][, ,] Value
            {
                get
                {
                    return this.value;
                }
                set
                {
                    this.value = value;
                }

            }

            ////////////////////
            /// Private Methods
            ////////////////////

            private Electric(int dimension, double[] window, int[] samplesize)
            {

                // Initialize field as zero everywhere

                int i, j, m, q;
                Dim = dimension;
                N = new int[dim];
                Size = new double[dim];
                N = samplesize;
                Size = window;
                for (i=0; i < dim; i++)
                    Value[i] = new double[N[0], N[1], N[2]];
                for (i = 0; i < N[0]; i++)
                {
                    for (j = 0; j < N[0]; j++)
                    {
                        for (m = 0; m < N[0]; m++)
                        {
                            for (q =0; q < dim; q++)
                                Value[q][i, j, m] = 0;
                        }
                    }
                }


            }



        }
        class Magnetic : Field<double[][,,]>
        {
            ////////////////////////
            //Property Declarations
            ///////////////////////


            private int dim;
            public int Dim
            {
                get
                {
                    return dim;
                }
                set
                {
                    dim = value;
                }

            }
            private double[] size;
            public double[] Size
            {
                get
                {
                    return size;
                }
                set
                {
                    size = value;
                }

            }
            private int[] n;
            public int[] N
            {
                get
                {
                    return n;
                }
                set
                {
                    n = value;
                }

            }
            private double[][,,] value;
            public double[][,,] Value
            {
                get
                {
                    return this.value;
                }
                set
                {
                    this.value = value;
                }

            }

            ////////////////////
            /// Private Methods
            ////////////////////

            private Magnetic(int dimension, double[] window, int[] samplesize)
            {

                // Initialize field as zero everywhere

                int i, j, m, q;
                Dim = dimension;
                N = new int[dim];
                Size = new double[dim];
                N = samplesize;
                Size = window;
                for (q = 0; q < dim; q++ )
                    Value[q] = new double[N[0], N[1], N[2]];
                for (i = 0; i < N[0]; i++ )
                {
                    for (j = 0; j < N[0]; j++)
                    {
                        for (m = 0; m < N[0]; m++)
                        {
                            for (q = 0; q < dim; q++)
                                Value[q][i, j, m] = 0;
                        }
                    }
                }
                    

            }



        }
        
        
}


