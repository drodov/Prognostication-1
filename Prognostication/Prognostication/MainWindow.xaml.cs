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
        Int32 K, N0 = 100;
        Double[] P;
        Double dt, M1, M2, S;
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
        }

        int CountNi(int idx)
        {
            idx--;
            if (idx == -1)
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
            if (i == -1)
                return 0;
            int i1 = i;
            int i2 = i + 1;
            return (double)(2.0 * ResCol[i].ni / (dt * (CountNi(i1) + CountNi(i2))));
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
        }

        void DrawErlGraphics()
        {
            double t, a1, a2, a3 = 0, sqr;
            a1 = 2.0 / M1;
            a2 = Math.Sqrt(6.0 / M2);
            for (int i = 0; i < K; i++)
            {
                sqr = Math.Sqrt((CountLambda(i) - CountLambda(i - 1)) / dt);
                a3 += (1.0 / K) * sqr / (1 - ((double)(i + 1) - 0.5) * dt * sqr);
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
        }

        void DrawVejbGraphics()
        {
            ZedGraphControl VejbZedGraph = VejbWFH.Child as ZedGraphControl;
            GraphPane pane = VejbZedGraph.GraphPane;
            pane.CurveList.Clear();
            PointPairList list = new PointPairList();
            // Заполняем список точек
            for (int i = 0; i < K; i++)
            {
                list.Add((ResCol[i].i - 0.5) * dt, P[i]);
            }
            LineItem myCurve = pane.AddCurve("Экспериментально", list, System.Drawing.Color.Red, SymbolType.None);
            VejbZedGraph.AxisChange();
            // Обновляем график
            VejbZedGraph.Invalidate();
        }

        void DrawNormGraphics()
        {
            double t, T1, T2, T3, T4, q1, q2, q3, q4;
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
            // Заполняем список точек
            for (int i = 0; i < K; i++)
            {
                t = (ResCol[i].i - 0.5) * dt;
                list0.Add(t, P[i]);
                list1.Add(t, FuncLowNorm(T1, q1, t));
                list2.Add(t, FuncLowNorm(T2, q2, t));
                list3.Add(t, FuncLowNorm(T3, q3, t));
            }
            LineItem Curve0 = pane.AddCurve("Экспериментально", list0, System.Drawing.Color.Red, SymbolType.None);
            LineItem Curve1 = pane.AddCurve("a51 = a52", list1, System.Drawing.Color.Black, SymbolType.None);
            //LineItem Curve2 = pane.AddCurve("a52", list2, System.Drawing.Color.Green, SymbolType.None);
            LineItem Curve3 = pane.AddCurve("a53", list3, System.Drawing.Color.Blue, SymbolType.None);
            NormZedGraph.AxisChange();
            // Обновляем график
            NormZedGraph.Invalidate();
        }

        void DrawShortNormGraphics()
        {
            ZedGraphControl ShortNormZedGraph = ShortNormWFH.Child as ZedGraphControl;
            GraphPane pane = ShortNormZedGraph.GraphPane;
            pane.CurveList.Clear(); 
            PointPairList list = new PointPairList();
            // Заполняем список точек
            for (int i = 0; i < K; i++)
            {
                list.Add((ResCol[i].i - 0.5) * dt, P[i]);
            }
            LineItem myCurve = pane.AddCurve("Экспериментально", list, System.Drawing.Color.Red, SymbolType.None);
            ShortNormZedGraph.AxisChange();
            // Обновляем график
            ShortNormZedGraph.Invalidate();
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
            return Math.Exp(-a * Math.Pow(t, 2));
        }

        double FuncLowVejb(double a, double b, double t)
        {
            return Math.Exp(-a * Math.Pow(t, b));
        }

        double FuncNorm(double T, double q, double t)
        {
            return (1 / (q * Math.Sqrt(2 * Math.PI))) * Math.Exp(-Math.Pow(t - T, 2) / (2 * Math.Pow(q, 2)));
        }

        double FuncShortNorm(double T, double q, double t)
        {
            return 0;
        }

        double FuncLowNorm(double T, double q, double t)
        {
            return 1 - IntegralForNorm(0, t, T, q);
        }

        double FuncLowShortNorm(double T, double q, double t)
        {
            return 0;
        }

        double IntegralForNorm(double lim1, double lim2, double T, double q)
        {
            double result = 0;
            double s = 0.0001;
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
    }
}
