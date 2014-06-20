using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using ILNumerics;
using ILNumerics.Drawing;
using ILNumerics.Drawing.Plotting;
using Fields;
using System.Threading;

namespace FDTDdemo
{
    public partial class Form1 : Form
    {
        Electric E;
        Magnetic H;
        int t = 0;
        int m;
        int[] samp = new int[3] { 200, 200, 200 };
        double[] wind = new double[3] { 200, 200, 200 };
        int imp0 = 377;
        ILArray<float> A = new ILArray<float>[200];
        ILScene scene = new ILScene();
        ILPlotCube plotCube;
       
        public Form1()
        {
            InitializeComponent();
            E = new Electric(3, wind, samp);
            H = new Magnetic(3, wind, samp);

            ilPanel1.Scene = scene;
            plotCube = scene.Add(new ILPlotCube());

            plotCube.Limits.Set(new Vector3(0f,-1f,0f), new Vector3(200f,1f,0f));
            plotCube.AllowRotation = false;
            plotCube.AllowZoom = false;
            plotCube.AllowPan = false;
            ilPanel1.Scene = scene;
            backgroundWorker1.RunWorkerAsync();
        }

        private void backgroundWorker1_DoWork(object sender, DoWorkEventArgs e)
        {
            int t;
            for (t = 0; t < 1000; t++)
            {
                //calculate fields
                

                for (m = 0; m < 199; m++)
                    H.Value[0][m, 0, 0] = H.Value[0][m, 0, 0] + (E.Value[0][m + 1, 0, 0] - E.Value[0][m, 0, 0]) / imp0;
                for (m = 1; m < 200; m++)
                    E.Value[0][m, 0, 0] = E.Value[0][m, 0, 0] + (H.Value[0][m, 0, 0] - H.Value[0][m - 1, 0, 0]) * imp0;
                
                E.Value[0][0, 0, 0] = Math.Exp(-(Math.Pow((t - 30f),  2f)) / 100f);
                

                for (m = 0; m < 200; m++)
                   A[m] = Convert.ToSingle(E.Value[0][m, 0, 0]);

            }
        }


        private void timer1_Tick(object sender, EventArgs e)
        {
            

                //calculate fields


                for (m = 0; m < 199; m++)
                    H.Value[0][m, 0, 0] = H.Value[0][m, 0, 0] + (E.Value[0][m + 1, 0, 0] - E.Value[0][m, 0, 0]) / imp0;
                for (m = 1; m < 200; m++)
                    E.Value[0][m, 0, 0] = E.Value[0][m, 0, 0] + (H.Value[0][m, 0, 0] - H.Value[0][m - 1, 0, 0]) * imp0;

                E.Value[0][0, 0, 0] = Math.Exp(-(Math.Pow((t - 30f), 2f)) / 100f);


                for (m = 0; m < 200; m++)
                    A[m] = Convert.ToSingle(E.Value[0][m, 0, 0]);

            
            plotCube.Children.Clear();
            plotCube.Add(new ILLinePlot(A.T));
            scene.Configure();
            t++;
        }

    }
}