using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FrameworkLibrary
{

        public interface Field
        {
            int Dim { get; set; }
            double[] Size { get; set; }
            int[] N { get; set; }
            var[][] Value { get; set; }
        }
        class Wavefunction : Field
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
        class Potential : Field
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

            private Potential(int dimension, double[] window, int[] samplesize)
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
}


