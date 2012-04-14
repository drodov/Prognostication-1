using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using System.Collections.ObjectModel;
using ZedGraph;

namespace Prognostication
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        double DELTA = 0.0001;
        Int32 K, N0 = 100;
        Double[] P;
        Double dt, M1, M2, S;
        Double[] D = new Double[16];
        ObservableCollection<Results> ResCol = new ObservableCollection<Results>();
        public MainWindow()
        {
            InitializeComponent();
            ZedGraphControl ExpZedGraph = ExpWFH.Child as ZedGraphControl;
            ZedGraphControl ErlZedGraph = ErlWFH.Child as ZedGraphControl;
            ZedGraphControl RelZedGraph = RelWFH.Child as ZedGraphControl;
            ZedGraphControl VejbZedGraph = VejbWFH.Child as ZedGraphControl;
            ZedGraphControl NormZedGraph = NormWFH.Child as ZedGraphControl;
            ZedGraphControl ShortNormZedGraph = ShortNormWFH.Child as ZedGraphControl;
            ExpZedGraph.GraphPane.Title.Text = "Экспоненциальный закон";
            ExpZedGraph.GraphPane.XAxis.Title.Text = "t";
            ExpZedGraph.GraphPane.YAxis.Title.Text = "P";
            ErlZedGraph.GraphPane.Title.Text = "Закон Эрланга";
            ErlZedGraph.GraphPane.XAxis.Title.Text = "t";
            ErlZedGraph.GraphPane.YAxis.Title.Text = "P";
            RelZedGraph.GraphPane.Title.Text = "Закон Рэлея";
            RelZedGraph.GraphPane.XAxis.Title.Text = "t";
            RelZedGraph.GraphPane.YAxis.Title.Text = "P";
            VejbZedGraph.GraphPane.Title.Text = "Закон Вейбулла";
            VejbZedGraph.GraphPane.XAxis.Title.Text = "t";
            VejbZedGraph.GraphPane.YAxis.Title.Text = "P";
            NormZedGraph.GraphPane.Title.Text = "Нормальный закон";
            NormZedGraph.GraphPane.XAxis.Title.Text = "t";
            NormZedGraph.GraphPane.YAxis.Title.Text = "P";
            ShortNormZedGraph.GraphPane.Title.Text = "Усеченный нормальный закон";
            ShortNormZedGraph.GraphPane.XAxis.Title.Text = "t";
            ShortNormZedGraph.GraphPane.YAxis.Title.Text = "P";
            dataGrid1.ItemsSource = ResCol;
        }

        private void KTextBox_TextChanged(object sender, TextChangedEventArgs e)
        {
            Int32.TryParse(KTextBox.Text, out K);
            if ( K < 10 )
                return;
            ResCol.Clear();
            for (int i = 0; i < K; i++)
            {
                ResCol.Add(new Results(i + 1));
            }
            dataGrid1.ItemsSource = ResCol;
        }

        private void dtTextBox_TextChanged(object sender, TextChangedEventArgs e)
        {
            Double.TryParse(dtTextBox.Text, out dt);
            if (K == 0)
                return;
        }

        private void StartButton_Click(object sender, RoutedEventArgs e)
        {
            M1 = 0;
            M2 = 0;
            S = 0;
            Double.TryParse(dtTextBox.Text, out dt);
            if (dt == 0)
                return;
            for (int i = 1; ; i *= 10 )
            {
                if ((int)dt / i == 0)
                {
                    DELTA *= i;
                    break;
                }
            }
            Int32.TryParse(KTextBox.Text, out K);
            if (K < 10)
                return;
            P = new Double[K];
            for(int i = 0; i < K; i++)
            {
                M1 += (double)(i + 1 - 0.5) * ResCol[i].ni;
                M2 += Math.Pow((double)(i + 1 - 0.5), 2) * ResCol[i].ni;
                P[i] = CountProbability(i);
            }
            M1 = M1 * dt / N0;
            M2 = M2 * Math.Pow( dt, 2 ) / N0;
            S = Math.Sqrt(M2 - Math.Pow(M1, 2));
            DrawExpGraphics();
            DrawErlGraphics();
            DrawRelGraphics();
            DrawVejbGraphics();
            DrawNormGraphics();
            DrawShortNormGraphics();
            ShowD();
        }

        int CountNi(int idx)
        {
            idx--;
            if (idx < 0)
                return N0;
            int Ni = N0;
            for (int i = 0; i < idx + 1; i++)
            {
                Ni -= ResCol[i].ni;
            }
            return Ni;
        }

        double CountProbability(int i)
        {
            int i1 = i;
            int i2 = i + 1;
            return (double)(CountNi(i1) + CountNi(i2)) / (2.0 * N0);
        }

        double CountLambda(int i)
        {
            if (i < 0)
                return 0;
            int i1 = i;
            int i2 = i + 1;
            return (double)(2.0 * ResCol[i].ni / (dt * (CountNi(i1) + CountNi(i2))));
        }

        double CountF(int i)
        {
            i--;
            return ((double)ResCol[i].ni / (N0 * dt));
        }

        void DrawExpGraphics()
        {
            double t, a1, a2, a3 = 0;
            a1 = 1.0 / M1;
            a2 = Math.Sqrt(2.0 / M2);
            for (int i = 0; i < K; i++)
            {
                a3 += ((double)(ResCol[i].ni)) / (CountNi(i) + CountNi(i + 1));
            }
            a3 = a3 * 2.0 / (K * dt);
            ZedGraphControl ExpZedGraph = ExpWFH.Child as ZedGraphControl;
            GraphPane pane = ExpZedGraph.GraphPane;
            pane.CurveList.Clear();
            PointPairList list0 = new PointPairList();
            PointPairList list1 = new PointPairList();
            PointPairList list2 = new PointPairList();
            PointPairList list3 = new PointPairList();
            // Заполняем список точек
            for (int i = 0; i < K; i++)
            {
                t = (ResCol[i].i - 0.5) * dt;
                list0.Add(t, P[i]);
                list1.Add(t, FuncLowExp(a1, t));
                list2.Add(t, FuncLowExp(a2, t));
                list3.Add(t, FuncLowExp(a3, t));
            }
            LineItem Curve0 = pane.AddCurve("Экспериментально", list0, System.Drawing.Color.Red, SymbolType.None);
            LineItem Curve1 = pane.AddCurve("a11", list1, System.Drawing.Color.Black, SymbolType.None);
            LineItem Curve2 = pane.AddCurve("a12", list2, System.Drawing.Color.Green, SymbolType.None);
            LineItem Curve3 = pane.AddCurve("a13", list3, System.Drawing.Color.Blue, SymbolType.None);
            ExpZedGraph.AxisChange();
            // Обновляем график
            ExpZedGraph.Invalidate();

            D[0] = CountDExp(a1);
            D[1] = CountDExp(a2);
            D[2] = CountDExp(a3);
        }

        void DrawErlGraphics()
        {
            double t, a1, a2, a3 = 0, sqr;
            a1 = 2.0 / M1;
            a2 = Math.Sqrt(6.0 / M2);
            for (int i = 0; i < K; i++)
            {
                sqr = Math.Sqrt((CountLambda(i) - CountLambda(i - 1)) / dt);/*
sqr = (CountLambda(i) - CountLambda(i - 1)) / dt;
if (sqr < 0)
    sqr *= -1;
sqr = Math.Sqrt(sqr)*/
                a3 += (1.0 / K) * sqr / (1 - (i + 1 - 0.5) * dt * sqr);
            }
            ZedGraphControl ErlZedGraph = ErlWFH.Child as ZedGraphControl;
            GraphPane pane = ErlZedGraph.GraphPane;
            pane.CurveList.Clear();
            PointPairList list0 = new PointPairList();
            PointPairList list1 = new PointPairList();
            PointPairList list2 = new PointPairList();
            PointPairList list3 = new PointPairList();
            // Заполняем список точек
            for (int i = 0; i < K; i++)
            {
                t = (ResCol[i].i - 0.5) * dt;
                list0.Add(t, P[i]);
                list1.Add(t, FuncLowErl(a1, t));
                list2.Add(t, FuncLowErl(a2, t));
                list3.Add(t, FuncLowErl(a3, t));
            }
            LineItem Curve0 = pane.AddCurve("Экспериментально", list0, System.Drawing.Color.Red, SymbolType.None);
            LineItem Curve1 = pane.AddCurve("a21", list1, System.Drawing.Color.Black, SymbolType.None);
            LineItem Curve2 = pane.AddCurve("a22", list2, System.Drawing.Color.Green, SymbolType.None);
            LineItem Curve3 = pane.AddCurve("a23", list3, System.Drawing.Color.Blue, SymbolType.None);
            ErlZedGraph.AxisChange();
            // Обновляем график
            ErlZedGraph.Invalidate();

            D[3] = CountDErl(a1);
            D[4] = CountDErl(a2);
            D[5] = CountDErl(a3);
        }

        void DrawRelGraphics()
        {
            double t, a1, a2, a3 = 0, a4;
            a1 = Math.PI / 4 / M1 / M1; 
            a2 = 1.0 / M2;
            a4 = 1.0 / (K * dt * dt) * (ResCol[K - 1].ni / CountNi(K - 1) - (double)ResCol[0].ni / (N0 + CountNi(1)));
            for (int i = 0; i < K; i++)
            {
                a3 += ((double)(ResCol[i].ni)) * (i + 1 - 0.5) / (CountNi(i) + CountNi(i + 1));
            }
            a3 = a3 * 12.0 / (K * dt * dt * (4 * K * K - 1));
            ZedGraphControl RelZedGraph = RelWFH.Child as ZedGraphControl;
            GraphPane pane = RelZedGraph.GraphPane;
            pane.CurveList.Clear();
            PointPairList list0 = new PointPairList();
            PointPairList list1 = new PointPairList();
            PointPairList list2 = new PointPairList();
            PointPairList list3 = new PointPairList();
            PointPairList list4 = new PointPairList();
            // Заполняем список точек
            for (int i = 0; i < K; i++)
            {
                t = (ResCol[i].i - 0.5) * dt;
                list0.Add(t, P[i]);
                list1.Add(t, FuncLowRel(a1, t));
                list2.Add(t, FuncLowRel(a2, t));
                list3.Add(t, FuncLowRel(a3, t));
                list4.Add(t, FuncLowRel(a4, t));
            }
            LineItem Curve0 = pane.AddCurve("Экспериментально", list0, System.Drawing.Color.Red, SymbolType.None);
            LineItem Curve1 = pane.AddCurve("a31", list1, System.Drawing.Color.Black, SymbolType.None);
            LineItem Curve2 = pane.AddCurve("a32", list2, System.Drawing.Color.Green, SymbolType.None);
            LineItem Curve3 = pane.AddCurve("a33", list3, System.Drawing.Color.Blue, SymbolType.None);
            LineItem Curve4 = pane.AddCurve("a34", list4, System.Drawing.Color.DarkOrange, SymbolType.None);
            RelZedGraph.AxisChange();
            // Обновляем график
            RelZedGraph.Invalidate();

            D[6] = CountDRel(a1);
            D[7] = CountDRel(a2);
            D[8] = CountDRel(a3);
            D[9] = CountDRel(a4);
        }

        void DrawVejbGraphics()
        {
            double t, x, y, a = 0, b = 0, c = 0, d = 0, e = 0;
            for (int i = 0; i < K; i++)
            {
                a += Math.Log(CountLambda(i)); // i+1
                b += Math.Log(i + 1 - 0.5) + Math.Log(dt);
                c += Math.Pow(Math.Log(i + 1 - 0.5) + Math.Log(dt), 2);
                e += Math.Log(CountLambda(i)) * (Math.Log(i + 1 - 0.5) + Math.Log(dt));
            }
            d = b * b;
            x = (a  * c - b * e) / (K * c - b * b);
            y = (a * b - K * e) / (b * b - K * c);
            double B = 1 + y;
            double A = Math.Exp(x) / B;
            ZedGraphControl VejbZedGraph = VejbWFH.Child as ZedGraphControl;
            GraphPane pane = VejbZedGraph.GraphPane;
            pane.CurveList.Clear();
            PointPairList list0 = new PointPairList();
            PointPairList list1 = new PointPairList();
            // Заполняем список точек
            for (int i = 0; i < K; i++)
            {
                t = (ResCol[i].i - 0.5) * dt;
                list0.Add(t, P[i]);
                list1.Add(t, FuncLowVejb(A, B, t));
            }//
            LineItem Curve0 = pane.AddCurve("Экспериментально", list0, System.Drawing.Color.Red, SymbolType.None);
            LineItem Curve1 = pane.AddCurve("a21", list1, System.Drawing.Color.Black, SymbolType.None);
            VejbZedGraph.AxisChange();
            // Обновляем график
            VejbZedGraph.Invalidate();

            D[10] = CountDVejb(A, B);
        }

        void DrawNormGraphics()
        {
            double t, T1, T2, T3, T4 = 0, q1, q2, q3, q4 = 0, a, b = 0, g = 1, Y, X;
            a = Math.Log(CountF(1) / CountF(K));
            for (int i = 1; i <= K - 1; i++)
            {
                b += Math.Pow(Math.Log(CountF(i) / CountF(i + 1)), 2);
                g *= Math.Pow(CountF(i) / CountF(i + 1), i);
            }
            g = Math.Log(g);
            Y = (double)(K - 1) * (0.5 * K * a - g) / (a * a - b * (K - 1));
            X = (0.5 * K * b * (K - 1) - g * a) / (b * (K - 1) - a * a);
            T4 = X * dt;
            q4 = Math.Sqrt(Y ) * dt;
            T1 = T2 = T3 = M1;
            q1 = q2 = S;
            q3 = q1 * Math.Sqrt(2);
            ZedGraphControl NormZedGraph = NormWFH.Child as ZedGraphControl;
            GraphPane pane = NormZedGraph.GraphPane;
            pane.CurveList.Clear();
            PointPairList list0 = new PointPairList();
            PointPairList list1 = new PointPairList();
            PointPairList list2 = new PointPairList();
            PointPairList list3 = new PointPairList();
            PointPairList list4 = new PointPairList();
            // Заполняем список точек
            for (int i = 0; i < K; i++)
            {
                t = (ResCol[i].i - 0.5) * dt;
                list0.Add(t, P[i]);
                list1.Add(t, FuncLowNorm(T1, q1, t));
                //                list2.Add(t, FuncLowNorm(T2, q2, t));
                list3.Add(t, FuncLowNorm(T3, q3, t));
                list4.Add(t, FuncLowNorm(T4, q4, t));
            }
            LineItem Curve0 = pane.AddCurve("Экспериментально", list0, System.Drawing.Color.Red, SymbolType.None);
            LineItem Curve1 = pane.AddCurve("a51 = a52", list1, System.Drawing.Color.Black, SymbolType.None);
            //LineItem Curve2 = pane.AddCurve("a52", list2, System.Drawing.Color.Green, SymbolType.None);
            LineItem Curve3 = pane.AddCurve("a53", list3, System.Drawing.Color.Blue, SymbolType.None);
            LineItem Curve4 = pane.AddCurve("a54", list4, System.Drawing.Color.Green, SymbolType.None);
            NormZedGraph.AxisChange();
            // Обновляем график
            NormZedGraph.Invalidate();

            D[11] = CountDNorm(q1, T1);
            D[12] = CountDNorm(q2, T2);
            D[13] = CountDNorm(q3, T3);
            D[14] = CountDNorm(q4, T4);
        }

        void DrawShortNormGraphics()
        {
            double t, C, T = 0, q = 0, a, b = 0, g = 1, Y, X;
            a = Math.Log(CountF(1) / CountF(K));
            for (int i = 1; i < K; i++)
            {
                b += Math.Pow(Math.Log(CountF(i) / CountF(i + 1)), 2);
                g *= Math.Pow(CountF(i) / CountF(i + 1), i);
            }
            g = Math.Log(g);
            Y = (double)(K - 1) * (0.5 * K * a - g) / (a * a - b * (K - 1));
            X = (0.5 * K * b * (K - 1) - g * a) / (b * (K - 1) - a * a);
            T = X * dt;
            q = Math.Sqrt(Y) * dt;
            C = 1.0 / (0.5 + FLaplas(T / q));
            ZedGraphControl ShortNormZedGraph = ShortNormWFH.Child as ZedGraphControl;
            GraphPane pane = ShortNormZedGraph.GraphPane;
            pane.CurveList.Clear();
            PointPairList list0 = new PointPairList();
            PointPairList list1 = new PointPairList();
            // Заполняем список точек
            for (int i = 0; i < K; i++)
            {
                t = (ResCol[i].i - 0.5) * dt;
                list0.Add(t, P[i]);
                list1.Add(t, FuncLowShortNorm(T, q, t, C));
            }
            LineItem Curve0 = pane.AddCurve("Экспериментально", list0, System.Drawing.Color.Red, SymbolType.None);
            LineItem Curve1 = pane.AddCurve("a6", list1, System.Drawing.Color.Black, SymbolType.None);
            ShortNormZedGraph.AxisChange();
            // Обновляем график
            ShortNormZedGraph.Invalidate();

            D[15] = CountDShortNorm(q, T, C);
        }

        double FuncLowExp(double a, double t)
        {
            return Math.Exp(-a * t);
        }

        double FuncLowErl(double a, double t)
        {
            return (1 + a * t) * Math.Exp(-a * t);
        }

        double FuncLowRel(double a, double t)
        {
            return Math.Exp(-a * t * t);
        }

        double FuncLowVejb(double a, double b, double t)
        {
            return Math.Exp(-a * Math.Pow(t, b));
        }

        double FuncNorm(double T, double q, double t)
        {
            return (1 / (q * Math.Sqrt(2 * Math.PI))) * Math.Exp(-Math.Pow(t - T, 2) / (2 * Math.Pow(q, 2)));
        }

        double FuncLowNorm(double T, double q, double t)
        {
            return 1 - IntegralForNorm(0, t, T, q);
        }

        double FuncLowShortNorm(double T, double q, double t, double C)
        {
            return C * (1 - IntegralForNorm(0, t, T, q));
        }

 /*       double FuncLowShortNorm(double T, double q, double t, double C)
        {
            double res = C * (0.5 - FLaplas((t - T) / q));
            return res;
        }
*/
        double IntegralForNorm(double lim1, double lim2, double T, double q)
        {
            double result = 0;
            double s = DELTA;
            for (double i = lim1 + s; i <= lim2 - s; i += 2 * s)
            {
                result += FuncNorm(T, q, i) * 2 * s;
            }
            return result;
        }

        double CountDExp(double a)
        {
            Double D = 0;
            Double temp;
            for (int i = 1; i <= K; i++)
            {
                temp = a * Math.Exp(-a * (i - 0.5) * dt) - ResCol[i - 1].ni / (N0 * dt);
                D += Math.Pow(temp, 2);
            }
            D /= K;
            return D;
        }

        double CountDErl(double a)
        {
            Double D = 0;
            Double temp;
            for (int i = 1; i <= K; i++)
            {
                temp = a * a * (i - 0.5) * dt * Math.Exp(-a * (i - 0.5) * dt) - ResCol[i - 1].ni / (N0 * dt);
                D += Math.Pow(temp, 2);
            }
            D /= K;
            return D;
        }

        double CountDRel(double a)
        {
            Double D = 0;
            Double temp;
            for (int i = 1; i <= K; i++)
            {
                temp = 2 * a * (i - 0.5) * dt * Math.Exp(-a * (i - 0.5) * (i - 0.5) * dt * dt) - ResCol[i - 1].ni / (N0 * dt);
                D += Math.Pow(temp, 2);
            }
            D /= K;
            return D;
        }

        double CountDVejb(double a, double b)
        {
            Double D = 0;
            Double temp;
            Double t;
            for (int i = 1; i <= K; i++)
            {
                t = (i - 0.5) * dt;
                temp = a * b * Math.Pow(t, b - 1) * Math.Exp(-a * Math.Pow(t, b)) - ResCol[i - 1].ni / (N0 * dt);
                D += Math.Pow(temp, 2);
            }
            D /= K;
            return D;
        }

        double CountDNorm(double q, double T)
        {
            Double D = 0;
            Double temp;
            Double t;
            for (int i = 1; i <= K; i++)
            {
                t = (i - 0.5) * dt;
                temp = 1 / (q * Math.Sqrt(2 * Math.PI)) * Math.Exp(-(t - T) * (t - T) / (2 * q * q)) - ResCol[i - 1].ni / (N0 * dt);
                D += Math.Pow(temp, 2);
            }
            D /= K;
            return D;
        }

        double CountDShortNorm(double q, double T, double C)
        {
            Double D = 0;
            Double temp;
            Double t;
            for (int i = 1; i <= K; i++)
            {
                t = (i - 0.5) * dt;
                temp = C / (q * Math.Sqrt(2 * Math.PI)) * Math.Exp(-(t - T) * (t - T) / (2 * q * q)) - ResCol[i - 1].ni / (N0 * dt);
                D += Math.Pow(temp, 2);
            }
            D /= K;
            return D;
        }

        void ShowD()
        {
            double min = D.Min();
            D11Label.Content += D[0].ToString();
            if (D[0] == min)
                D11Label.Foreground = Brushes.Red;
            D12Label.Content += D[1].ToString();
            if (D[1] == min)
                D12Label.Foreground = Brushes.Red;
            D13Label.Content += D[2].ToString();
            if (D[2] == min)
                D13Label.Foreground = Brushes.Red;
            D21Label.Content += D[3].ToString();
            if (D[3] == min)
                D21Label.Foreground = Brushes.Red;
            D22Label.Content += D[4].ToString();
            if (D[4] == min)
                D22Label.Foreground = Brushes.Red;
            D23Label.Content += D[5].ToString();
            if (D[5] == min)
                D23Label.Foreground = Brushes.Red;
            D31Label.Content += D[6].ToString();
            if (D[6] == min)
                D31Label.Foreground = Brushes.Red;
            D32Label.Content += D[7].ToString();
            if (D[7] == min)
                D32Label.Foreground = Brushes.Red;
            D33Label.Content += D[8].ToString();
            if (D[8] == min)
                D33Label.Foreground = Brushes.Red;
            D34Label.Content += D[9].ToString();
            if (D[9] == min)
                D34Label.Foreground = Brushes.Red;
            D4Label.Content += D[10].ToString();
            if (D[10] == min)
                D4Label.Foreground = Brushes.Red;
            D51Label.Content += D[11].ToString();
            if (D[11] == min)
                D51Label.Foreground = Brushes.Red;
            D52Label.Content += D[12].ToString();
            if (D[12] == min)
                D52Label.Foreground = Brushes.Red;
            D53Label.Content += D[13].ToString();
            if (D[13] == min)
                D53Label.Foreground = Brushes.Red;
            D54Label.Content += D[14].ToString();
            if (D[14] == min)
                D54Label.Foreground = Brushes.Red;
            D6Label.Content += D[15].ToString();
            if (D[15] == min)
                D6Label.Foreground = Brushes.Red;
            AnswDLabel.Content = min.ToString();
        }

        double Factorial(int x)
        {
            double y = 1;
            if (x == 0)
                return 0;
            for (int i = x; i != 0; i--)
            {
                y *= i;
            }
                return y;
        }

        double FLaplas(double x)
        {
            double res = 1.0 / Math.Sqrt(2 * Math.PI) * IntegralForLaplas(0, x, x);
            return res;
        }

        double IntegralForLaplas(double lim1, double lim2, double x)
        {/*
            if (lim2 < 0)
            {
                double temp = lim2;
                lim2 = lim1;
                lim1 = temp;
            }*/
            double result = 0;
            double s = DELTA;
            for (double i = lim1 + s; i <= lim2 - s; i += 2 * s)
            {
                result += Math.Exp(-(x * x) / 2) * s * 2;
            }
            return result;
        }
    }
}
