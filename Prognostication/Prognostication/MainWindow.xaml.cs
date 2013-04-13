using System;
using System.Linq;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Media;
using System.Collections.ObjectModel;
using ZedGraph;

namespace Prognostication
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        private double _delta;
        private int K;
        private const int N0 = 100;
        private double[] P;
        private double dt;
        private double M1;
        private double M2;
        private double S;
        private double[] D = new double[22];  // D*
        private double[] DD = new double[22]; // D
        private ObservableCollection<Results> _resultCollection = new ObservableCollection<Results>();
        private readonly ObservableCollection<Results> _tempResultCollection = new ObservableCollection<Results>();

        public MainWindow()
        {
            InitializeComponent();

            var expZedGraph = (ZedGraphControl)ExpWFH.Child;
            expZedGraph.GraphPane.Title.Text = "Экспоненциальный закон";
            expZedGraph.GraphPane.XAxis.Title.Text = "t";
            expZedGraph.GraphPane.YAxis.Title.Text = "P";
            
            var erlZedGraph = (ZedGraphControl)ErlWFH.Child;
            erlZedGraph.GraphPane.Title.Text = "Закон Эрланга";
            erlZedGraph.GraphPane.XAxis.Title.Text = "t";
            erlZedGraph.GraphPane.YAxis.Title.Text = "P";
            
            var relZedGraph = (ZedGraphControl)RelWFH.Child;
            relZedGraph.GraphPane.Title.Text = "Закон Рэлея";
            relZedGraph.GraphPane.XAxis.Title.Text = "t";
            relZedGraph.GraphPane.YAxis.Title.Text = "P";

            var vejbZedGraph = (ZedGraphControl)VejbWFH.Child;
            vejbZedGraph.GraphPane.Title.Text = "Закон Вейбулла";
            vejbZedGraph.GraphPane.XAxis.Title.Text = "t";
            vejbZedGraph.GraphPane.YAxis.Title.Text = "P";

            var normZedGraph = (ZedGraphControl)NormWFH.Child;
            normZedGraph.GraphPane.Title.Text = "Нормальный закон";
            normZedGraph.GraphPane.XAxis.Title.Text = "t";
            normZedGraph.GraphPane.YAxis.Title.Text = "P";

            var shortNormZedGraph = (ZedGraphControl)ShortNormWFH.Child;
            shortNormZedGraph.GraphPane.Title.Text = "Усеченный нормальный закон";
            shortNormZedGraph.GraphPane.XAxis.Title.Text = "t";
            shortNormZedGraph.GraphPane.YAxis.Title.Text = "P";

            var logNormZedGraph = (ZedGraphControl)LogNormWFH.Child;
            logNormZedGraph.GraphPane.Title.Text = "Логнормальный закон";
            logNormZedGraph.GraphPane.XAxis.Title.Text = "t";
            logNormZedGraph.GraphPane.YAxis.Title.Text = "P";

            var gammaZedGraph = (ZedGraphControl)GammaWFH.Child;
            gammaZedGraph.GraphPane.Title.Text = "Гамма закон";
            gammaZedGraph.GraphPane.XAxis.Title.Text = "t";
            gammaZedGraph.GraphPane.YAxis.Title.Text = "P";
            
            dataGrid1.ItemsSource = _resultCollection;
        }

        private void KTextBox_TextChanged(object sender, TextChangedEventArgs e)
        {
            Int32.TryParse(KTextBox.Text, out K);
            if (K < 1)
            {
                return;
            }

            // save old items
            for (int i = 0; i < _resultCollection.Count; i++ )
            {
                if(i < _tempResultCollection.Count)
                {
                    _tempResultCollection[i] = _resultCollection[i];
                }
                else
                {
                    _tempResultCollection.Add(_resultCollection[i]);
                }
            }

            // new items
            _resultCollection = new ObservableCollection<Results>();
            for (int i = 1; i <= K; i++)
            {
                _resultCollection.Add(_tempResultCollection.Count >= i
                               ? new Results(i, _tempResultCollection[i - 1].Ni)
                               : new Results(i));
            }
            dataGrid1.ItemsSource = _resultCollection;
        }

        private void StartButton_Click(object sender, RoutedEventArgs e)
        {
            try
            {
                Clear();
                M1 = 0;
                M2 = 0;
                S = 0;
                double ni = _resultCollection.Sum(res => res.Ni);
                if (ni != N0)
                {
                    MessageBox.Show(string.Format("Сумма ni ({0}) не равна N0 (100).", ni));
                    return;
                }
                double.TryParse(dtTextBox.Text, out dt);
                if (dt == 0 || dt <= 0)
                {
                    MessageBox.Show("Неправильно указано время dt.\n(Пример: 10,1)");
                    return;
                }
                _delta = 0.0001;
                for (int i = 1; ; i *= 10)
                {
                    if ((int)dt / i == 0)
                    {
                        _delta *= i;
                        break;
                    }
                }
                Int32.TryParse(KTextBox.Text, out K);
                if (K < 5)
                {
                    MessageBox.Show("Неправильно указано число интервалов K либо оно мало.");
                    return;
                }
                P = new double[K];
                for (int i = 0; i < K; i++)
                {
                    M1 += (i + 1 - 0.5) * _resultCollection[i].Ni;
                    M2 += Math.Pow(i + 1 - 0.5, 2) * _resultCollection[i].Ni;
                    P[i] = CountProbability(i);
                }
                M1 = M1 * dt / N0;
                M2 = M2 * Math.Pow(dt, 2) / N0;
                S = Math.Sqrt(M2 - Math.Pow(M1, 2));
                DrawExpGraphics();
                DrawErlGraphics();
                DrawRelGraphics();
                DrawVejbGraphics();
                DrawNormGraphics();
                DrawShortNormGraphics();
                DrawLogNormGraphics();
                DrawGammaGraphics();
                tabControl1.SelectedIndex = 9;
                ShowD();
            }
            catch (Exception ex)
            {
                MessageBox.Show(ex.Message);
                Clear();
            }
        }

        private double CountNi(int idx)
        {
            idx--;
            if (idx < 0)
                return N0;
            double Ni = N0;
            for (int i = 0; i < idx + 1; i++)
            {
                Ni -= _resultCollection[i].Ni;
            }
            return Ni;
        }

        private double CountProbability(int i)
        {
            int i1 = i;
            int i2 = i + 1;
            return (CountNi(i1) + CountNi(i2)) / (2.0 * N0);
        }

        private double CountLambda(int i)
        {
            if (i < 0)
            {
                return 0;
            }
            int i1 = i;
            int i2 = i + 1;
            return 2.0 * _resultCollection[i].Ni / (dt * (CountNi(i1) + CountNi(i2)));
        }

        private double CountF(int i)
        {
            i--;
            return _resultCollection[i].Ni / (N0 * dt);
        }

        #region Draw graphics

        private void DrawExpGraphics()
        {
            double a3 = 0;
            double a1 = 1.0 / M1;
            double a2 = Math.Sqrt(2.0 / M2);
            for (int i = 0; i < K; i++)
            {
                a3 += _resultCollection[i].Ni / (CountNi(i) + CountNi(i + 1));
            }
            a3 = a3 * 2.0 / (K * dt);
            var expZedGraph = (ZedGraphControl)ExpWFH.Child;
            GraphPane pane = expZedGraph.GraphPane;
            pane.CurveList.Clear();
            var list0 = new PointPairList();
            var list1 = new PointPairList();
            var list2 = new PointPairList();
            var list3 = new PointPairList();
            // Заполняем список точек
            for (int i = 0; i < K; i++)
            {
                double t = (_resultCollection[i].I - 0.5) * dt;
                list0.Add(t, P[i]);
                list1.Add(t, FuncLowExp(a1, t));
                list2.Add(t, FuncLowExp(a2, t));
                list3.Add(t, FuncLowExp(a3, t));
            }
            pane.AddCurve("Экспериментально", list0, System.Drawing.Color.Red, SymbolType.None);
            pane.AddCurve("a11", list1, System.Drawing.Color.Black, SymbolType.None);
            pane.AddCurve("a12", list2, System.Drawing.Color.Green, SymbolType.None);
            pane.AddCurve("a13", list3, System.Drawing.Color.Blue, SymbolType.None);
            expZedGraph.AxisChange();
            // Обновляем график
            expZedGraph.Invalidate();

            if (double.IsInfinity(a1) || double.IsNaN(a1))
            {
                A11Label.Content = "a11 = --";
                D[0] = 100;
                DD[0] = 100;
            }
            else
            {
                A11Label.Content = "a11 = " + a1.ToString();
                D[0] = CalculateDExp(a1);
                DD[0] = CalculateDDExp(a1);
            }
            if (double.IsInfinity(a2) || double.IsNaN(a2))
            {
                A12Label.Content = "a12 = --";
                D[1] = 100;
                DD[1] = 100;
            }
            else
            {
                A12Label.Content = "a12 = " + a2.ToString();
                D[1] = CalculateDExp(a2);
                DD[1] = CalculateDDExp(a2);
            }
            if (double.IsInfinity(a3) || double.IsNaN(a3))
            {
                A13Label.Content = "a13 = --";
                D[2] = 100;
                DD[2] = 100;
            }
            else
            {
                A13Label.Content = "a13 = " + a3.ToString();
                D[2] = CalculateDExp(a3);
                DD[2] = CalculateDDExp(a3);
            }
        }

        private void DrawErlGraphics()
        {
            double a1 = 2.0 / M1;
            
            double a2 = Math.Sqrt(6.0 / M2);

            double a3 = 0;
            for (int i = 0; i < K; i++)
            {
                double sqr = Math.Sqrt((CountLambda(i) - CountLambda(i - 1)) / dt);
                a3 += sqr / Math.Abs((1 - (i + 1 - 0.5) * dt * sqr));
            }
            a3 /= K;
// статья
            //double a4 = 0;
            //for (int i = 0; i < K - 1; i++)
            //{
            //    a4 += (CountLambda(i) / (i + 1 - 0.5) - (CountLambda(i + 1) - CountLambda(i))) / ((i + 1 + 0.5) * (CountLambda(i + 1) - CountLambda(i)));
            //}
            //a4 /= ((K - 1) * dt);
//
            var erlZedGraph = (ZedGraphControl)ErlWFH.Child;
            GraphPane pane = erlZedGraph.GraphPane;
            pane.CurveList.Clear();
            var list0 = new PointPairList();
            var list1 = new PointPairList();
            var list2 = new PointPairList();
            var list3 = new PointPairList();
//статья
//            var list4 = new PointPairList();
            // Заполняем список точек
            for (int i = 0; i < K; i++)
            {
                double t = (_resultCollection[i].I - 0.5) * dt;
                list0.Add(t, P[i]);
                list1.Add(t, FuncLowErl(a1, t));
                list2.Add(t, FuncLowErl(a2, t));
                list3.Add(t, FuncLowErl(a3, t));
//статья
//                list4.Add(t, FuncLowErl(a4, t));
            }
            pane.AddCurve("Экспериментально", list0, System.Drawing.Color.Red, SymbolType.None);
            pane.AddCurve("a21", list1, System.Drawing.Color.Black, SymbolType.None);
            pane.AddCurve("a22", list2, System.Drawing.Color.Green, SymbolType.None);
            pane.AddCurve("a23", list3, System.Drawing.Color.Blue, SymbolType.None);
//статья
//            pane.AddCurve("a24", list4, System.Drawing.Color.Blue, SymbolType.None);
            erlZedGraph.AxisChange();
            // Обновляем график
            erlZedGraph.Invalidate();

            if (double.IsInfinity(a1) || double.IsNaN(a1))
            {
                A21Label.Content = "a21 = --";
                D[3] = 100;
                DD[3] = 100;
            }
            else
            {
                A21Label.Content = "a21 = " + a1.ToString();
                D[3] = CalculateDErl(a1);
                DD[3] = CalculateDDErl(a1);
            }
            if (double.IsInfinity(a2) || double.IsNaN(a2))
            {
                A22Label.Content = "a22 = --";
                D[4] = 100;
                DD[4] = 100;
            }
            else
            {
                A22Label.Content = "a22 = " + a2.ToString();
                D[4] = CalculateDErl(a2);
                DD[4] = CalculateDDErl(a2);
            }
            if (double.IsInfinity(a3) || double.IsNaN(a3))
            {
                A23Label.Content = "a23 = --";
                D[5] = 100;
                DD[5] = 100;
            }
            else
            {
                A23Label.Content = "a23 = " + a3.ToString();
                D[5] = CalculateDErl(a3);
                DD[5] = CalculateDDErl(a3);
            }
        }

        private void DrawRelGraphics()
        {
            double a3 = 0;
            double a1 = Math.PI / 4 / M1 / M1;
            double a2 = 1.0 / M2;
            double a4 = 1.0 / (K * dt * dt) * (_resultCollection[K - 1].Ni / CountNi(K - 1) - _resultCollection[0].Ni / (N0 + CountNi(1)));
            for (int i = 0; i < K; i++)
            {
                a3 += _resultCollection[i].Ni * (i + 1 - 0.5) / (CountNi(i) + CountNi(i + 1));
            }
            a3 = a3 * 12.0 / (K * dt * dt * (4 * K * K - 1));
            var relZedGraph = (ZedGraphControl)RelWFH.Child;
            GraphPane pane = relZedGraph.GraphPane;
            pane.CurveList.Clear();
            var list0 = new PointPairList();
            var list1 = new PointPairList();
            var list2 = new PointPairList();
            var list3 = new PointPairList();
            var list4 = new PointPairList();
            // Заполняем список точек
            for (int i = 0; i < K; i++)
            {
                double t = (_resultCollection[i].I - 0.5) * dt;
                list0.Add(t, P[i]);
                list1.Add(t, FuncLowRel(a1, t));
                list2.Add(t, FuncLowRel(a2, t));
                list3.Add(t, FuncLowRel(a3, t));
                list4.Add(t, FuncLowRel(a4, t));
            }
            pane.AddCurve("Экспериментально", list0, System.Drawing.Color.Red, SymbolType.None);
            pane.AddCurve("a31", list1, System.Drawing.Color.Black, SymbolType.None);
            pane.AddCurve("a32", list2, System.Drawing.Color.Green, SymbolType.None);
            pane.AddCurve("a33", list3, System.Drawing.Color.Blue, SymbolType.None);
            pane.AddCurve("a34", list4, System.Drawing.Color.DarkOrange, SymbolType.None);
            relZedGraph.AxisChange();
            // Обновляем график
            relZedGraph.Invalidate();

            if (double.IsInfinity(a1) || double.IsNaN(a1))
            {
                A31Label.Content = "a31 = --";
                D[6] = 100;
                DD[6] = 100;
            }
            else
            {
                A31Label.Content = "a31 = " + a1.ToString();
                D[6] = CalculateDRel(a3);
                DD[6] = CalculateDDRel(a3);
            }
            if (double.IsInfinity(a2) || double.IsNaN(a2))
            {
                A32Label.Content = "a32 = --";
                D[7] = 100;
                DD[7] = 100;
            }
            else
            {
                A32Label.Content = "a32 = " + a2.ToString();
                D[7] = CalculateDRel(a2);
                DD[7] = CalculateDDRel(a2);
            }
            if (double.IsInfinity(a3) || double.IsNaN(a3))
            {
                A33Label.Content = "a33 = --";
                D[8] = 100;
                DD[8] = 100;
            }
            else
            {
                A33Label.Content = "a33 = " + a3.ToString();
                D[8] = CalculateDRel(a3);
                DD[8] = CalculateDDRel(a3);
            }
            if (double.IsInfinity(a4) || double.IsNaN(a4))
            {
                A34Label.Content = "a34 = --";
                D[9] = 100;
                DD[9] = 100;
            }
            else
            {
                A34Label.Content = "a34 = " + a4.ToString();
                D[9] = CalculateDRel(a4);
                DD[9] = CalculateDDRel(a4);
            }
        }

        private void DrawVejbGraphics()
        {
            double a = 0, b = 0, c = 0, e = 0;
            for (int i = 0; i < K; i++)
            {
                a += Math.Log(CountLambda(i)); // i+1
                b += Math.Log(i + 1 - 0.5) + Math.Log(dt);
                c += Math.Pow(Math.Log(i + 1 - 0.5) + Math.Log(dt), 2);
                e += Math.Log(CountLambda(i)) * (Math.Log(i + 1 - 0.5) + Math.Log(dt));
            }
            double x = (a * c - b * e) / (K * c - b * b);
            double y = (a * b - K * e) / (b * b - K * c);
            double B = 1 + y;
            double A = Math.Exp(x) / B;
            var VejbZedGraph = (ZedGraphControl)VejbWFH.Child;
            GraphPane pane = VejbZedGraph.GraphPane;
            pane.CurveList.Clear();
            var list0 = new PointPairList();
            var list1 = new PointPairList();
            // Заполняем список точек
            for (int i = 0; i < K; i++)
            {
                double t = (_resultCollection[i].I - 0.5) * dt;
                list0.Add(t, P[i]);
                list1.Add(t, FuncLowVejb(A, B, t));
            }
            //
            pane.AddCurve("Экспериментально", list0, System.Drawing.Color.Red, SymbolType.None);
            pane.AddCurve("a21", list1, System.Drawing.Color.Black, SymbolType.None);
            VejbZedGraph.AxisChange();
            // Обновляем график
            VejbZedGraph.Invalidate();

            if (double.IsNaN(A) || double.IsInfinity(A) || double.IsNaN(B) || double.IsInfinity(B))
            {
                ALabel.Content = "a = --";
                BLabel.Content = "b = --";
                D[10] = 100;
                DD[10] = 100;
            }
            else
            {
                ALabel.Content = "a = " + A.ToString();
                BLabel.Content = "b = " + B.ToString();
                D[10] = CalculateDVejb(A, B);
                DD[10] = CalculateDDVejb(A, B);
            }
        }

        private void DrawNormGraphics()
        {
            double T2, T3, q2, b = 0, g = 1;
            double a = Math.Log(CountF(1) / CountF(K));
            for (int i = 1; i <= K - 1; i++)
            {
                b += Math.Pow(Math.Log(CountF(i) / CountF(i + 1)), 2);
                g *= Math.Pow(CountF(i) / CountF(i + 1), i);
            }
            g = Math.Log(g);
            double Y = (K - 1) * (0.5 * K * a - g) / (a * a - b * (K - 1));
            double X = (0.5 * K * b * (K - 1) - g * a) / (b * (K - 1) - a * a);
            double T4 = X * dt;
            double q4 = Math.Sqrt(Y) * dt;
            double T1 = T2 = T3 = M1;
            double q1 = q2 = S;
            double q3 = q1 * Math.Sqrt(2);
            var normZedGraph = (ZedGraphControl)NormWFH.Child;
            GraphPane pane = normZedGraph.GraphPane;
            pane.CurveList.Clear();
            var list0 = new PointPairList();
            var list1 = new PointPairList();
            // var list2 = new PointPairList();
            var list3 = new PointPairList();
            var list4 = new PointPairList();
            // Заполняем список точек
            for (int i = 0; i < K; i++)
            {
                double t = (_resultCollection[i].I - 0.5) * dt;
                list0.Add(t, P[i]);
                list1.Add(t, FuncLowNorm(T1, q1, t));
                //                list2.Add(t, FuncLowNorm(T2, q2, t));
                list3.Add(t, FuncLowNorm(T3, q3, t));
                list4.Add(t, FuncLowNorm(T4, q4, t));
            }
            pane.AddCurve("Экспериментально", list0, System.Drawing.Color.Red, SymbolType.None);
            pane.AddCurve("a51 = a52", list1, System.Drawing.Color.Black, SymbolType.None);
            //pane.AddCurve("a52", list2, System.Drawing.Color.Green, SymbolType.None);
            pane.AddCurve("a53", list3, System.Drawing.Color.Blue, SymbolType.None);
            pane.AddCurve("a54", list4, System.Drawing.Color.Green, SymbolType.None);
            normZedGraph.AxisChange();
            // Обновляем график
            normZedGraph.Invalidate();

            if (double.IsNaN(q1) || double.IsInfinity(q1) || double.IsNaN(T1) || double.IsInfinity(T1))
            {
                T52Label.Content = "T = --";
                Q52Label.Content = "q = --";
                D[11] = 100;
                D[12] = 100;
                DD[11] = 100;
                DD[12] = 100;
            }
            else
            {
                T52Label.Content = "T = " + T2.ToString();
                Q52Label.Content = "q = " + q2.ToString();
                D[11] = CalculateDNorm(q1, T1);
                D[12] = CalculateDNorm(q2, T2);
                DD[11] = CalculateDDNorm(q1, T1);
                DD[12] = CalculateDDNorm(q2, T2);
            }
            if (double.IsNaN(q3) || double.IsInfinity(q3) || double.IsNaN(T3) || double.IsInfinity(T3))
            {
                T53Label.Content = "T = --";
                Q53Label.Content = "q = --";
                D[13] = 100;
                DD[13] = 100;
            }
            else
            {
                T53Label.Content = "T = " + T3.ToString();
                Q53Label.Content = "q = " + q3.ToString();
                D[13] = CalculateDNorm(q3, T3);
                DD[13] = CalculateDDNorm(q3, T3);
            }
            if (double.IsNaN(q4) || double.IsInfinity(q4) || double.IsNaN(T4) || double.IsInfinity(T4))
            {
                T54Label.Content = "T = --";
                Q54Label.Content = "q = --";
                D[14] = 100;
                DD[14] = 100;
            }
            else
            {
                T54Label.Content = "T = " + T4.ToString();
                Q54Label.Content = "q = " + q4.ToString();
                D[14] = CalculateDNorm(q4, T4);
                DD[14] = CalculateDDNorm(q4, T4);
            }
        }

        private void DrawShortNormGraphics()
        {
            double b = 0, g = 1;
            double a = Math.Log(CountF(1) / CountF(K));
            for (int i = 1; i < K; i++)
            {
                b += Math.Pow(Math.Log(CountF(i) / CountF(i + 1)), 2);
                g *= Math.Pow(CountF(i) / CountF(i + 1), i);
            }
            g = Math.Log(g);
            double Y = (K - 1) * (0.5 * K * a - g) / (a * a - b * (K - 1));
            double X = (0.5 * K * b * (K - 1) - g * a) / (b * (K - 1) - a * a);
            double T = X * dt;
            double q = Math.Sqrt(Y) * dt;
            double C = 1.0 / (0.5 + FLaplas(T / q));
            var shortNormZedGraph = (ZedGraphControl)ShortNormWFH.Child;
            GraphPane pane = shortNormZedGraph.GraphPane;
            pane.CurveList.Clear();
            var list0 = new PointPairList();
            var list1 = new PointPairList();
            // Заполняем список точек
            for (int i = 0; i < K; i++)
            {
                double t = (_resultCollection[i].I - 0.5) * dt;
                list0.Add(t, P[i]);
                list1.Add(t, FuncLowShortNorm(T, q, t, C));
            }
            pane.AddCurve("Экспериментально", list0, System.Drawing.Color.Red, SymbolType.None);
            pane.AddCurve("a6", list1, System.Drawing.Color.Black, SymbolType.None);
            shortNormZedGraph.AxisChange();
            // Обновляем график
            shortNormZedGraph.Invalidate();

            if (double.IsNaN(q) || double.IsInfinity(q) || double.IsNaN(T) || double.IsInfinity(T))
            {
                T6Label.Content = "T = --";
                Q6Label.Content = "q = --";
                D[15] = 100;
                DD[15] = 100;
            }
            else
            {
                T6Label.Content = "T = " + T.ToString();
                Q6Label.Content = "q = " + q.ToString();
                D[15] = CalculateDShortNorm(q, T, C);
                DD[15] = CalculateDDShortNorm(q, T, C);
            }
        }

        private void DrawLogNormGraphics()
        {
            // 1
            double q1 = Math.Sqrt(Math.Log(M2 / M1 / M1, Math.E));
            double m1 = 0.5 * Math.Log(M1 * M1 * M1 * M1 / M2, Math.E);

            // 2
            double q2 = 0;
            double m2 = 0;
            for (int i = 0; i < K; i++)
            {
                double lnTi = Math.Log((_resultCollection[i].I - 0.5) * dt, Math.E);
                m2 += lnTi;
                q2 += lnTi * lnTi;
            }
            m2 /= K;
            q2 /= K;
            q2 -= m2;

            // 3
            double q = 0, p = 0, z = 0;
            for (int i = 1; i <= K - 1; i++)
            {
                double num = Math.Log(CountF(i + 1) / CountF(i)) / Math.Log((i + 0.5) / (i - 0.5));
                p += num;
                q += num * num;
                z += num * Math.Log(i * i - 0.25) / 2;
            }
            p += K - 1;
            q += 2 * p + K - 1;
            double y = (K - 1) * Math.Log(dt / 4) + Math.Log(Factorial(2 * K - 2) * Math.Sqrt(2 * K - 1) / Factorial(K - 1));
            z += y + (p + 1 - K) * Math.Log(dt);
            double m3 = (y * q - p * z) / ((K - 1) * q - p * p);
            double q3 = (y * p - (K - 1) * z) / ((K - 1) * q - p * p);

            // draw
            var logNormZedGraph = (ZedGraphControl)LogNormWFH.Child;
            GraphPane pane = logNormZedGraph.GraphPane;
            pane.CurveList.Clear();
            var list0 = new PointPairList();
            var list1 = new PointPairList();
            var list2 = new PointPairList();
            var list3 = new PointPairList();
            // Заполняем список точек
            for (int i = 0; i < K; i++)
            {
                double t = (_resultCollection[i].I - 0.5) * dt;
                list0.Add(t, P[i]);
                list1.Add(t, FuncLowLogNorm(t, m1, q1));
                list2.Add(t, FuncLowLogNorm(t, m2, q2));
                list3.Add(t, FuncLowLogNorm(t, m3, q3));
            }
            pane.AddCurve("Экспериментально", list0, System.Drawing.Color.Red, SymbolType.None);
            pane.AddCurve("a1", list1, System.Drawing.Color.Black, SymbolType.None);
            pane.AddCurve("a2", list2, System.Drawing.Color.Blue, SymbolType.None);
            pane.AddCurve("a3", list3, System.Drawing.Color.Green, SymbolType.None);
            logNormZedGraph.AxisChange();
            // Обновляем график
            logNormZedGraph.Invalidate();

            if (double.IsNaN(q1) || double.IsInfinity(q1) || q1 == 0 || double.IsNaN(m1) || double.IsInfinity(m1))
            {
                Q71Label.Content = "q = --";
                M71Label.Content = "m = --";
                D[16] = 100;
                DD[16] = 100;
            }
            else
            {
                Q71Label.Content = "q = " + q1.ToString();
                M71Label.Content = "m = " + m1.ToString();
                D[16] = CalculateDLogNorm(m1, q1);
                DD[16] = CalculateDDLogNorm(m1, q1);
            }
            if (double.IsNaN(q2) || double.IsInfinity(q2) || q2 == 0 || double.IsNaN(m2) || double.IsInfinity(m2))
            {
                Q72Label.Content = "q = --";
                M72Label.Content = "m = --";
                D[17] = 100;
                DD[17] = 100;
            }
            else
            {
                Q72Label.Content = "q = " + q2.ToString();
                M72Label.Content = "m = " + m2.ToString();
                D[17] = CalculateDLogNorm(m2, q2);
                DD[17] = CalculateDDLogNorm(m2, q2);
            }
            if (double.IsNaN(q3) || double.IsInfinity(q3) || q3 == 0 || double.IsNaN(m3) || double.IsInfinity(m3))
            {
                Q73Label.Content = "q = --";
                M73Label.Content = "m = --";
                D[18] = 100;
                DD[18] = 100;
            }
            else
            {
                Q73Label.Content = "q = " + q3.ToString();
                M73Label.Content = "m = " + m3.ToString();
                D[18] = CalculateDLogNorm(m3, q3);
                DD[18] = CalculateDDLogNorm(m3, q3);
            }
        }

        private void DrawGammaGraphics()
        {
            // 1
            double a1 = M1*M1/(M2 - M1*M1);
            double b1 = (M2 - M1*M1)/M1;

            // 2
            const double a2 = 2.25;
            double b2 = K*dt/4.5;

            // 3
            double num1 = 0, num2 = 0;
            for (int i = 1; i <= K - 1; i++)
            {
                double tempLogValue = Math.Log((i + 0.5)/(i - 0.5));
                num1 += Math.Log(CountF(i + 1) / CountF(i)) * tempLogValue;
                num2 += tempLogValue*tempLogValue;
            }
            double tempCalculationValue = (Math.Log(2*K - 1)*Math.Log(CountF(K)/CountF(1)) - (K - 1)*num1)/
                                          (Math.Pow(Math.Log(2*K - 1), 2) - (K - 1)*num2);
            double a3 = 1 + tempCalculationValue;
            double b3 = dt/tempCalculationValue;

            // draw
            var gammaZedGraph = (ZedGraphControl)GammaWFH.Child;
            GraphPane pane = gammaZedGraph.GraphPane;
            pane.CurveList.Clear();
            var list0 = new PointPairList();
            var list1 = new PointPairList();
            var list2 = new PointPairList();
            var list3 = new PointPairList();
            // Заполняем список точек
            for (int i = 0; i < K; i++)
            {
                double t = (_resultCollection[i].I - 0.5) * dt;
                list0.Add(t, P[i]);
                list1.Add(t, FuncLowGamma(t, a1, b1));
                list2.Add(t, FuncLowGamma(t, a2, b2));
                list3.Add(t, FuncLowGamma(t, a3, b3));
            }
            pane.AddCurve("Экспериментально", list0, System.Drawing.Color.Red, SymbolType.None);
            pane.AddCurve("a1", list1, System.Drawing.Color.Black, SymbolType.None);
            pane.AddCurve("a2", list2, System.Drawing.Color.Blue, SymbolType.None);
            pane.AddCurve("a3", list3, System.Drawing.Color.Green, SymbolType.None);
            gammaZedGraph.AxisChange();
            // Обновляем график
            gammaZedGraph.Invalidate();

            if (double.IsNaN(a1) || double.IsInfinity(a1) || double.IsNaN(b1) || double.IsInfinity(b1))
            {
                A81Label.Content = "a = --";
                B81Label.Content = "b = --";
                D[19] = 100;
                DD[19] = 100;
            }
            else
            {
                A81Label.Content = "a = " + a1.ToString();
                B81Label.Content = "b = " + b1.ToString();
                D[19] = CalculateDGamma(a1, b1);
                DD[19] = CalculateDDGamma(a1, b1);
            }
            if (double.IsNaN(a2) || double.IsInfinity(a2) || double.IsNaN(b2) || double.IsInfinity(b2))
            {
                A82Label.Content = "a = --";
                B82Label.Content = "b = --";
                D[20] = 100;
                DD[20] = 100;
            }
            else
            {
                A82Label.Content = "a = " + a2.ToString();
                B82Label.Content = "b = " + b2.ToString();
                D[20] = CalculateDGamma(a2, b2);
                DD[20] = CalculateDDGamma(a2, b2);
            }
            if (double.IsNaN(a3) || double.IsInfinity(a3) || double.IsNaN(b3) || double.IsInfinity(b3))
            {
                A83Label.Content = "a = --";
                B83Label.Content = "b = --";
                D[21] = 100;
                DD[21] = 100;
            }
            else
            {
                A83Label.Content = "a = " + a3.ToString();
                B83Label.Content = "b = " + b3.ToString();
                D[21] = CalculateDGamma(a3, b3);
                DD[21] = CalculateDDGamma(a3, b3);
            }
        }

        #endregion

        #region Lows

        private double FuncLowExp(double a, double t)
        {
            return Math.Exp(-a * t);
        }

        private double FuncLowErl(double a, double t)
        {
            return (1 + a * t) * Math.Exp(-a * t);
        }

        private double FuncLowRel(double a, double t)
        {
            return Math.Exp(-a * t * t);
        }

        private double FuncLowVejb(double a, double b, double t)
        {
            return Math.Exp(-a * Math.Pow(t, b));
        }

        private double FuncLowNorm(double T, double q, double t)
        {
            return 1 - IntegralForNorm(0, t, T, q);
        }

        private double FuncLowShortNorm(double T, double q, double t, double C)
        {
            return C * (1 - IntegralForNorm(0, t, T, q));
        }

        /*       double FuncLowShortNorm(double T, double q, double t, double C)
               {
                   double res = C * (0.5 - FLaplas((t - T) / q));
                   return res;
               }
       */

        private double FuncLowLogNorm(double t, double m, double q)
        {
            return q == 0 ? double.NaN : (0.5 - FLaplas((Math.Log(t, Math.E) - m)/q));
        }

        private double IntegralForNorm(double lim1, double lim2, double T, double q)
        {
            double result = 0;
            double s = _delta;
            for (double i = lim1 + s; i <= lim2 - s; i += 2 * s)
            {
                result += FuncNorm(T, q, i) * 2 * s;
            }
            return result;
        }

        private double FuncNorm(double T, double q, double t)
        {
            return (1 / (q * Math.Sqrt(2 * Math.PI))) * Math.Exp(-Math.Pow(t - T, 2) / (2 * Math.Pow(q, 2)));
        }

        private double FuncLowGamma(double t, double a, double b)
        {
            return 1 - IntegralForGamma(0, t, a, b);
        }

        private double IntegralForGamma(double lim1, double lim2, double a, double b)
        {
            double result = 0;
            double s = _delta;
            for (double i = lim1 + s; i <= lim2 - s; i += 2 * s)
            {
                result += FuncGamma(a, b, i) * 2 * s;
            }
            return result;
        }

        private double FuncGamma(double a, double b, double t)
        {
            return 1/(Math.Pow(b, a)*Gamma(a))*Math.Pow(t, a - 1)*Math.Exp(-t/b);
        }

        private readonly double[] _cof =
            {
                2.5066282746310005,
                1.0000000000190015,
                76.18009172947146,
                -86.50532032941677,
                24.01409824083091,
                -1.231739572450155,
                0.1208650973866179e-2,
                -0.5395239384953e-5
            };

        /* logarithm of gamma-function */
        private double GammLn(double x)
        {
            /* calculate the series */
            double ser = _cof[1];
            double y = x;
            for (int j = 2; j < 8; j++)
            {
                y += 1.0;
                ser += _cof[j]/y;
            }
            /* and the other parts of the function */
            y = x + 5.5;
            y -= (x + 0.5)*Math.Log(y);
            return (-y + Math.Log(_cof[0]*ser/x));
        }

        /* the gamma-function itself */
        private double Gamma(double x)
        {
            return Math.Exp(GammLn(x));
        }

        #endregion

        #region Calculate D*

        private double CalculateDExp(double a)
        {
            double D = 0;
            for (int i = 1; i <= K; i++)
            {
                double temp = a * Math.Exp(-a * (i - 0.5) * dt) - _resultCollection[i - 1].Ni / (N0 * dt);
                D += Math.Pow(temp, 2);
            }
            D /= K;
            if (double.IsNaN(D) || double.IsInfinity(D))
            {
                D = 100;
            }
            return D;
            /*
            double D = 0;
            double temp;
            for (int i = 1; i <= K; i++)
            {
                temp = Math.Exp(-a * (i - 0.5) * dt) - 0.5 * (CountProbability(i - 2) + CountProbability(i - 1));
                D += Math.Pow(temp, 2);
            }
            D /= K;
            if (D.ToString() == "NaN" || D.ToString() == "Infinity")
                D = 100;
            return D;*/
        }

        private double CalculateDErl(double a)
        {
            double D = 0;
            for (int i = 1; i <= K; i++)
            {
                double temp = a * a * (i - 0.5) * dt * Math.Exp(-a * (i - 0.5) * dt) - _resultCollection[i - 1].Ni / (N0 * dt);
                D += Math.Pow(temp, 2);
            }
            D /= K;
            if (double.IsNaN(D) || double.IsInfinity(D))
            {
                D = 100;
            }
            return D;
            /*
            double D = 0;
            double temp;
            for (int i = 1; i <= K; i++)
            {
                temp = (1 + a * (i - 0.5) * dt) * Math.Exp(-a * (i - 0.5) * dt) - 0.5 * (CountProbability(i - 2) + CountProbability(i - 1));
                D += Math.Pow(temp, 2);
            }
            D /= K;
            if (D.ToString() == "NaN" || D.ToString() == "Infinity")
                D = 100;
            return D;*/
        }

        private double CalculateDRel(double a)
        {
            double D = 0;
            for (int i = 1; i <= K; i++)
            {
                double temp = 2 * a * (i - 0.5) * dt * Math.Exp(-a * (i - 0.5) * (i - 0.5) * dt * dt) - _resultCollection[i - 1].Ni / (N0 * dt);
                D += Math.Pow(temp, 2);
            }
            D /= K;
            if (double.IsNaN(D) || double.IsInfinity(D))
            {
                D = 100;
            }
            return D;
            /*
            double D = 0;
            double temp;
            for (int i = 1; i <= K; i++)
            {
                temp = Math.Exp(-a * (i - 0.5) * (i - 0.5) * dt * dt) - 0.5 * (CountProbability(i - 2) + CountProbability(i - 1));
                D += Math.Pow(temp, 2);
            }
            D /= K;
            if (D.ToString() == "NaN" || D.ToString() == "Infinity")
                D = 100;
            return D;*/
        }

        private double CalculateDVejb(double a, double b)
        {
            double D = 0;
            for (int i = 1; i <= K; i++)
            {
                double t = (i - 0.5) * dt;
                double temp = a * b * Math.Pow(t, b - 1) * Math.Exp(-a * Math.Pow(t, b)) - _resultCollection[i - 1].Ni / (N0 * dt);
                D += Math.Pow(temp, 2);
            }
            D /= K;
            if (double.IsNaN(D) || double.IsInfinity(D))
            {
                D = 100;
            }
            return D;
            /*
            double D = 0;
            double temp;
            for (int i = 1; i <= K; i++)
            {
                temp = Math.Exp(-a * Math.Pow(((i - 0.5) * dt), b)) - 0.5 * (CountProbability(i - 2) + CountProbability(i - 1));
                D += Math.Pow(temp, 2);
            }
            D /= K;
            if (D.ToString() == "NaN" || D.ToString() == "Infinity")
                D = 100;
            return D;*/
        }

        private double CalculateDNorm(double q, double T)
        {
            double D = 0;
            for (int i = 1; i <= K; i++)
            {
                double t = (i - 0.5) * dt;
                double temp = 1 / (q * Math.Sqrt(2 * Math.PI)) * Math.Exp(-(t - T) * (t - T) / (2 * q * q)) - _resultCollection[i - 1].Ni / (N0 * dt);
                D += Math.Pow(temp, 2);
            }
            D /= K;
            if (double.IsNaN(D) || double.IsInfinity(D))
            {
                D = 100;
            }
            return D;
            /*
            double D = 0;
            double temp;
            for (int i = 1; i <= K; i++)
            {
                temp = 0.5 - FLaplas(((i - 0.5) * dt - T) / q) - 0.5 * (CountProbability(i - 2) + CountProbability(i - 1));
                D += Math.Pow(temp, 2);
            }
            D /= K;
            if (D.ToString() == "NaN" || D.ToString() == "Infinity")
                D = 100;
            return D;*/
        }

        private double CalculateDShortNorm(double q, double T, double C)
        {
            double D = 0;
            for (int i = 1; i <= K; i++)
            {
                double t = (i - 0.5) * dt;
                double temp = C/(q*Math.Sqrt(2*Math.PI))*Math.Exp(-(t - T)*(t - T)/(2*q*q)) -
                              _resultCollection[i - 1].Ni/(N0*dt);
                D += Math.Pow(temp, 2);
            }
            D /= K;
            if (double.IsNaN(D) || double.IsInfinity(D))
            {
                D = 100;
            }
            return D;
            /*
            double D = 0;
            double temp;
            for (int i = 1; i <= K; i++)
            {
                temp = C * (0.5 - FLaplas(((i - 0.5) * dt - T) / q)) - 0.5 * (CountProbability(i - 2) + CountProbability(i - 1));
                D += Math.Pow(temp, 2);
            }
            D /= K;
            if (D.ToString() == "NaN" || D.ToString() == "Infinity")
                D = 100;
            return D;*/
        }

        private double CalculateDLogNorm(double m, double q)
        {
            double d = 0;
            for (int i = 1; i <= K; i++)
            {
                double t = (i - 0.5) * dt;
                d += Math.Pow(CountProbability(i - 1) - FuncLowLogNorm(t, m, q), 2);
            }
            d /= K;
            if (double.IsNaN(d) || double.IsInfinity(d))
            {
                d = 100;
            }
            return d;
        }

        private double CalculateDGamma(double a, double b)
        {
            double d = 0;
            for (int i = 1; i <= K; i++)
            {
                double t = (i - 0.5) * dt;
                d += Math.Pow(CountProbability(i - 1) - FuncLowGamma(t, a, b), 2);
            }
            d /= K;
            if (double.IsNaN(d) || double.IsInfinity(d))
            {
                d = 100;
            }
            return d;
        }

        #endregion

        #region Calculate D

        private double CalculateDDExp(double a)
        {
            double D = 0;
            for (int i = 1; i <= K; i++)
            {
                double temp = Math.Exp(-a * (i - 0.5) * dt) - 0.5 * (CountProbability(i - 2) + CountProbability(i - 1));
                D += Math.Pow(temp, 2);
            }
            D /= K;
            if (double.IsNaN(D) || double.IsInfinity(D))
            {
                D = 100;
            }
            return D;
        }

        private double CalculateDDErl(double a)
        {
            double D = 0;
            for (int i = 1; i <= K; i++)
            {
                double temp = (1 + a * (i - 0.5) * dt) * Math.Exp(-a * (i - 0.5) * dt) - 0.5 * (CountProbability(i - 2) + CountProbability(i - 1));
                D += Math.Pow(temp, 2);
            }
            D /= K;
            if (double.IsNaN(D) || double.IsInfinity(D))
            {
                D = 100;
            }
            return D;
        }

        private double CalculateDDRel(double a)
        {
            double D = 0;
            for (int i = 1; i <= K; i++)
            {
                double temp = Math.Exp(-a * (i - 0.5) * (i - 0.5) * dt * dt) - 0.5 * (CountProbability(i - 2) + CountProbability(i - 1));
                D += Math.Pow(temp, 2);
            }
            D /= K;
            if (double.IsNaN(D) || double.IsInfinity(D))
            {
                D = 100;
            }
            return D;
        }

        private double CalculateDDVejb(double a, double b)
        {
            double D = 0;
            for (int i = 1; i <= K; i++)
            {
                double temp = Math.Exp(-a * Math.Pow(((i - 0.5) * dt), b)) - 0.5 * (CountProbability(i - 2) + CountProbability(i - 1));
                D += Math.Pow(temp, 2);
            }
            D /= K;
            if (double.IsNaN(D) || double.IsInfinity(D))
            {
                D = 100;
            }
            return D;
        }

        private double CalculateDDNorm(double q, double T)
        {
            double D = 0;
            for (int i = 1; i <= K; i++)
            {
                double temp = 0.5 - FLaplas(((i - 0.5) * dt - T) / q) - 0.5 * (CountProbability(i - 2) + CountProbability(i - 1));
                D += Math.Pow(temp, 2);
            }
            D /= K;
            if (double.IsNaN(D) || double.IsInfinity(D))
            {
                D = 100;
            }
            return D;
        }

        private double CalculateDDShortNorm(double q, double T, double C)
        {
            double D = 0;
            for (int i = 1; i <= K; i++)
            {
                double temp = C * (0.5 - FLaplas(((i - 0.5) * dt - T) / q)) - 0.5 * (CountProbability(i - 2) + CountProbability(i - 1));
                D += Math.Pow(temp, 2);
            }
            D /= K;
            if (double.IsInfinity(D) || double.IsNaN(D))
            {
                D = 100;
            }
            return D;
        }

        private double CalculateDDLogNorm(double m, double q)
        {
            double d = 0;
            for (int i = 1; i <= K; i++)
            {
                double t = (i - 0.5) * dt;
                d += Math.Pow(CountF(i) - FuncLowLogNorm(t, m, q), 2);
            }
            d /= K;
            if (double.IsNaN(d) || double.IsInfinity(d))
            {
                d = 100;
            }
            return d;
        }

        private double CalculateDDGamma(double a, double b)
        {
            double d = 0;
            for (int i = 1; i <= K; i++)
            {
                double t = (i - 0.5) * dt;
                d += Math.Pow(CountF(i) - FuncLowGamma(t, a, b), 2);
            }
            d /= K;
            if (double.IsNaN(d) || double.IsInfinity(d))
            {
                d = 100;
            }
            return d;
        }

        #endregion

        private void ShowD()
        {
            string ans = string.Empty;
            double min = D.Min();
            bool flagAnalys = false;
            bool flagChangLow = false;

            #region Exponent

            if (D[0] == 100)
            {
                D11Label.Content += "--";
                DD11Label.Content += "--";
                A11Label.Content = string.Empty;
            }
            else
            {
                D11Label.Content += D[0].ToString();
                DD11Label.Content += DD[0].ToString();
            }
            if (D[0] == min)
            {
                ans = "Экспоненциальное";
                D11Label.Foreground = Brushes.Red;
            }

            if (D[1] == 100)
            {
                D12Label.Content += "--";
                DD12Label.Content += "--";
                A12Label.Content = string.Empty;
            }
            else
            {
                D12Label.Content += D[1].ToString();
                DD12Label.Content += DD[1].ToString();
            }
            if (D[1] == min)
            {
                ans = "Экспоненциальное";
                D12Label.Foreground = Brushes.Red;
            }

            if (D[2] == 100)
            {
                D13Label.Content += "--";
                DD13Label.Content += "--";
                A13Label.Content = string.Empty;
            }
            else
            {
                D13Label.Content += D[2].ToString();
                DD13Label.Content += DD[2].ToString();
            }
            if (D[2] == min)
            {
                ans = "Экспоненциальное";
                D13Label.Foreground = Brushes.Red;
            }

            #endregion
            #region Erlang

            if (D[3] == 100)
            {
                D21Label.Content += "--";
                DD21Label.Content += "--";
                A21Label.Content = string.Empty;
            }
            else
            {
                D21Label.Content += D[3].ToString();
                DD21Label.Content += DD[3].ToString();
            }
            if (D[3] == min)
            {
                ans = "Эрланга";
                D21Label.Foreground = Brushes.Red;
            }

            if (D[4] == 100)
            {
                D22Label.Content += "--";
                DD22Label.Content += "--";
                A22Label.Content = string.Empty;
            }
            else
            {
                D22Label.Content += D[4].ToString();
                DD22Label.Content += DD[4].ToString();
            }
            if (D[4] == min)
            {
                ans = "Эрланга";
                D22Label.Foreground = Brushes.Red;
            }

            if (D[5] == 100)
            {
                D23Label.Content += "--";
                DD23Label.Content += "--";
                A23Label.Content = string.Empty;
            }
            else
            {
                D23Label.Content += D[5].ToString();
                DD23Label.Content += DD[5].ToString();
            }
            if (D[5] == min)
            {
                ans = "Эрланга";
                D23Label.Foreground = Brushes.Red;
            }

            #endregion
            #region Relej

            if (D[6] == 100)
            {
                D31Label.Content += "--";
                DD31Label.Content += "--";
                A31Label.Content = string.Empty;
            }
            else
            {
                D31Label.Content += D[6].ToString();
                DD31Label.Content += DD[6].ToString();
            }
            if (D[6] == min)
            {
                ans = "Рэлея";
                D31Label.Foreground = Brushes.Red;
            }

            if (D[7] == 100)
            {
                D32Label.Content += "--";
                DD32Label.Content += "--";
                A32Label.Content = string.Empty;
            }
            else
            {
                D32Label.Content += D[7].ToString();
                DD32Label.Content += DD[7].ToString();
            }
            if (D[7] == min)
            {
                ans = "Рэлея";
                D32Label.Foreground = Brushes.Red;
            }

            if (D[8] == 100)
            {
                D33Label.Content += "--";
                DD33Label.Content += "--";
                A33Label.Content = string.Empty;
            }
            else
            {
                D33Label.Content += D[8].ToString();
                DD33Label.Content += DD[8].ToString();
            }
            if (D[8] == min)
            {
                ans = "Рэлея";
                D33Label.Foreground = Brushes.Red;
            }

            if (D[9] == 100)
            {
                D34Label.Content += "--";
                DD34Label.Content += "--";
                A34Label.Content = string.Empty;
            }
            else
            {
                D34Label.Content += D[9].ToString();
                DD34Label.Content += DD[9].ToString();
            }
            if (D[9] == min)
            {
                ans = "Рэлея";
                D34Label.Foreground = Brushes.Red;
            }

            #endregion
            #region Vejbul

            if (D[10] == 100)
            {
                D4Label.Content += "--";
                DD4Label.Content += "--";
                ALabel.Content = string.Empty;
                BLabel.Content = string.Empty;
            }
            else
            {
                D4Label.Content += D[10].ToString();
                DD4Label.Content += DD[10].ToString();
            }
            if (D[10] == min)
            {
                ans = "Вейбулла";
                D4Label.Foreground = Brushes.Red;
            }

            #endregion
            #region Normal

            if (D[11] == 100)
            {
                D51Label.Content += "--";
                DD51Label.Content += "--";
            }
            else
            {
                D51Label.Content += D[11].ToString();
                DD51Label.Content += DD[11].ToString();
            }
            if (D[11] == min)
            {
                ans = "Нормальный";
                D51Label.Foreground = Brushes.Red;
                flagAnalys = true;
            }

            if (D[12] == 100)
            {
                D52Label.Content += "--";
                DD52Label.Content += "--";
                T52Label.Content = string.Empty;
                Q52Label.Content = string.Empty;
            }
            else
            {
                D52Label.Content += D[12].ToString();
                DD52Label.Content += DD[12].ToString();
            }
            if (D[12] == min)
            {
                ans = "Нормальный";
                D52Label.Foreground = Brushes.Red;
                flagAnalys = true;
            }

            if (D[13] == 100)
            {
                D53Label.Content += "--";
                DD53Label.Content += "--";
                T53Label.Content = string.Empty;
                Q53Label.Content = string.Empty;
            }
            else
            {
                D53Label.Content += D[13].ToString();
                DD53Label.Content += DD[13].ToString();
            }
            if (D[13] == min)
            {
                ans = "Нормальный";
                D53Label.Foreground = Brushes.Red;
                flagAnalys = true;
            }

            if (D[14] == 100)
            {
                D54Label.Content += "--";
                DD54Label.Content += "--";
                T54Label.Content = string.Empty;
                Q54Label.Content = string.Empty;
            }
            else
            {
                D54Label.Content += D[14].ToString();
                DD54Label.Content += DD[14].ToString();
            }
            if (D[14] == min)
            {
                ans = "Нормальный";
                D54Label.Foreground = Brushes.Red;
                flagAnalys = true;
            }

            #endregion
            #region Short Norm

            if (D[15] == 100)
            {
                D6Label.Content += "--";
                DD6Label.Content += "--";
                T6Label.Content = string.Empty;
                Q6Label.Content = string.Empty;
            }
            else
            {
                D6Label.Content += D[15].ToString();
                DD6Label.Content += DD[15].ToString();
            }
            if (D[15] == min)
            {
                ans = "Усеч. нормальный";
                D6Label.Foreground = Brushes.Red;
            }

            #endregion
            #region Log Norm

            if (D[16] == 100)
            {
                D71Label.Content += "--";
                DD71Label.Content += "--";
                Q71Label.Content = string.Empty;
                M71Label.Content = string.Empty;
            }
            else
            {
                D71Label.Content += D[16].ToString();
                DD71Label.Content += DD[16].ToString();
            }
            if (D[16] == min)
            {
                ans = "Логнормальный";
                D71Label.Foreground = Brushes.Red;
            }

            if (D[17] == 100)
            {
                D72Label.Content += "--";
                DD72Label.Content += "--";
                Q72Label.Content = string.Empty;
                M72Label.Content = string.Empty;
            }
            else
            {
                D72Label.Content += D[17].ToString();
                DD72Label.Content += DD[17].ToString();
            }
            if (D[17] == min)
            {
                ans = "Логнормальный";
                D72Label.Foreground = Brushes.Red;
            }

            if (D[18] == 100)
            {
                D73Label.Content += "--";
                DD73Label.Content += "--";
                Q73Label.Content = string.Empty;
                M73Label.Content = string.Empty;
            }
            else
            {
                D73Label.Content += D[18].ToString();
                DD73Label.Content += DD[18].ToString();
            }
            if (D[18] == min)
            {
                ans = "Логнормальный";
                D73Label.Foreground = Brushes.Red;
            }

            #endregion
            #region Gamma

            if (D[19] == 100)
            {
                D81Label.Content += "--";
                DD81Label.Content += "--";
                A81Label.Content = string.Empty;
                B81Label.Content = string.Empty;
            }
            else
            {
                D81Label.Content += D[19].ToString();
                DD81Label.Content += DD[19].ToString();
            }
            if (D[19] == min)
            {
                ans = "Гамма";
                D81Label.Foreground = Brushes.Red;
            }

            if (D[20] == 100)
            {
                D82Label.Content += "--";
                DD82Label.Content += "--";
                A82Label.Content = string.Empty;
                B82Label.Content = string.Empty;
            }
            else
            {
                D82Label.Content += D[20].ToString();
                DD82Label.Content += DD[20].ToString();
            }
            if (D[20] == min)
            {
                ans = "Гамма";
                D82Label.Foreground = Brushes.Red;
            }

            if (D[21] == 100)
            {
                D83Label.Content += "--";
                DD83Label.Content += "--";
                A83Label.Content = string.Empty;
                B83Label.Content = string.Empty;
            }
            else
            {
                D83Label.Content += D[21].ToString();
                DD83Label.Content += DD[21].ToString();
            }
            if (D[21] == min)
            {
                ans = "Гамма";
                D83Label.Foreground = Brushes.Red;
            }

            #endregion

            if (flagAnalys)
            {
                flagChangLow = Analys();
            }

            if (flagChangLow)
            {
                ans = "Усеч. нормальный";
                MessageBox.Show("По D - нормальный закон.\nОднако, судя по данным Ni - это усеченный нормальный.");
            }

            AnswDLabel.Content = string.Format("D* = {0}", min);
            AnswLowLabel.Content = ans;
        }

        private double FLaplas(double x)
        {
            return 1.0 / Math.Sqrt(2 * Math.PI) * IntegralForLaplas(0, x, x);
        }

        private double IntegralForLaplas(double lim1, double lim2, double x)
        {
            /*
            if (lim2 < 0)
            {
                double temp = lim2;
                lim2 = lim1;
                lim1 = temp;
            }*/
            double result = 0;
            double s = _delta;
            for (double i = lim1 + s; i <= lim2 - s; i += 2 * s)
            {
                result += Math.Exp(-(x * x) / 2) * s * 2;
            }
            return result;
        }
            
        private double Factorial(int x)
        {
            double y = 1;
            if (x == 0)
            {
                return 0;
            }

            for (int i = x; i != 0; i--)
            {
                y *= i;
            }
            return y;
        }

        private void Clear()
        {
            var expZedGraph = (ZedGraphControl)ExpWFH.Child;
            var erlZedGraph = (ZedGraphControl)ErlWFH.Child;
            var relZedGraph = (ZedGraphControl)RelWFH.Child;
            var vejbZedGraph = (ZedGraphControl)VejbWFH.Child;
            var normZedGraph = (ZedGraphControl)NormWFH.Child;
            var shortNormZedGraph = (ZedGraphControl)ShortNormWFH.Child;
            var logNormZedGraph = (ZedGraphControl)NormWFH.Child;
            var gammaZedGraph = (ZedGraphControl)ShortNormWFH.Child;
            
            var pane = expZedGraph.GraphPane;
            pane.CurveList.Clear();
            pane = erlZedGraph.GraphPane;
            pane.CurveList.Clear();
            pane = relZedGraph.GraphPane;
            pane.CurveList.Clear();
            pane = vejbZedGraph.GraphPane;
            pane.CurveList.Clear();
            pane = normZedGraph.GraphPane;
            pane.CurveList.Clear();
            pane = shortNormZedGraph.GraphPane;
            pane.CurveList.Clear();
            pane = logNormZedGraph.GraphPane;
            pane.CurveList.Clear();
            pane = gammaZedGraph.GraphPane;
            pane.CurveList.Clear();


            D11Label.Content = "D*11 = ";
            D11Label.Foreground = Brushes.Black;
            D12Label.Content = "D*12 = ";
            D12Label.Foreground = Brushes.Black;
            D13Label.Content = "D*13 = ";
            D13Label.Foreground = Brushes.Black;
            D21Label.Content = "D*21 = ";
            D21Label.Foreground = Brushes.Black;
            D22Label.Content = "D*22 = ";
            D22Label.Foreground = Brushes.Black;
            D23Label.Content = "D*23 = ";
            D23Label.Foreground = Brushes.Black;
            D31Label.Content = "D*31 = ";
            D31Label.Foreground = Brushes.Black;
            D32Label.Content = "D*32 = ";
            D32Label.Foreground = Brushes.Black;
            D33Label.Content = "D*33 = ";
            D33Label.Foreground = Brushes.Black;
            D34Label.Content = "D*34 = ";
            D34Label.Foreground = Brushes.Black;
            D4Label.Content = "D*4 = ";
            D4Label.Foreground = Brushes.Black;
            D51Label.Content = "D*51 = ";
            D51Label.Foreground = Brushes.Black;
            D52Label.Content = "D*52 = ";
            D52Label.Foreground = Brushes.Black;
            D53Label.Content = "D*53 = ";
            D53Label.Foreground = Brushes.Black;
            D54Label.Content = "D*54 = ";
            D54Label.Foreground = Brushes.Black;
            D6Label.Content = "D*6 = ";
            D6Label.Foreground = Brushes.Black;
            D71Label.Content = "D*71 = ";
            D71Label.Foreground = Brushes.Black;
            D72Label.Content = "D*72 = ";
            D72Label.Foreground = Brushes.Black;
            D73Label.Content = "D*73 = ";
            D73Label.Foreground = Brushes.Black;
            D81Label.Content = "D*81 = ";
            D81Label.Foreground = Brushes.Black;
            D82Label.Content = "D*82 = ";
            D82Label.Foreground = Brushes.Black;
            D83Label.Content = "D*83 = ";
            D83Label.Foreground = Brushes.Black;

            DD11Label.Content = "D11 = ";
            DD11Label.Foreground = Brushes.Black;
            DD12Label.Content = "D12 = ";
            DD12Label.Foreground = Brushes.Black;
            DD13Label.Content = "D13 = ";
            DD13Label.Foreground = Brushes.Black;
            DD21Label.Content = "D21 = ";
            DD21Label.Foreground = Brushes.Black;
            DD22Label.Content = "D22 = ";
            DD22Label.Foreground = Brushes.Black;
            DD23Label.Content = "D23 = ";
            DD23Label.Foreground = Brushes.Black;
            DD31Label.Content = "D31 = ";
            DD31Label.Foreground = Brushes.Black;
            DD32Label.Content = "D32 = ";
            DD32Label.Foreground = Brushes.Black;
            DD33Label.Content = "D33 = ";
            DD33Label.Foreground = Brushes.Black;
            DD34Label.Content = "D34 = ";
            DD34Label.Foreground = Brushes.Black;
            DD4Label.Content = "D4 = ";
            DD4Label.Foreground = Brushes.Black;
            DD51Label.Content = "D51 = ";
            DD51Label.Foreground = Brushes.Black;
            DD52Label.Content = "D52 = ";
            DD52Label.Foreground = Brushes.Black;
            DD53Label.Content = "D53 = ";
            DD53Label.Foreground = Brushes.Black;
            DD54Label.Content = "D54 = ";
            DD54Label.Foreground = Brushes.Black;
            DD6Label.Content = "D6 = ";
            DD6Label.Foreground = Brushes.Black;
            DD71Label.Content = "D71 = ";
            DD71Label.Foreground = Brushes.Black;
            DD72Label.Content = "D72 = ";
            DD72Label.Foreground = Brushes.Black;
            DD73Label.Content = "D73 = ";
            DD73Label.Foreground = Brushes.Black;
            DD81Label.Content = "D81 = ";
            DD81Label.Foreground = Brushes.Black;
            DD82Label.Content = "D82 = ";
            DD82Label.Foreground = Brushes.Black;
            DD83Label.Content = "D83 = ";
            DD83Label.Foreground = Brushes.Black;
            
            AnswDLabel.Content = string.Empty;
            AnswLowLabel.Content = string.Empty;

            A11Label.Content = string.Empty;
            A12Label.Content = string.Empty;
            A13Label.Content = string.Empty;
            
            A21Label.Content = string.Empty;
            A22Label.Content = string.Empty;
            A23Label.Content = string.Empty;
            
            A31Label.Content = string.Empty;
            A32Label.Content = string.Empty;
            A33Label.Content = string.Empty;
            A34Label.Content = string.Empty;
            
            T52Label.Content = string.Empty;
            Q52Label.Content = string.Empty;
            T53Label.Content = string.Empty;
            Q53Label.Content = string.Empty;
            T54Label.Content = string.Empty;
            Q54Label.Content = string.Empty;
            
            T6Label.Content = string.Empty;
            Q6Label.Content = string.Empty;
            
            Q71Label.Content = string.Empty;
            M71Label.Content = string.Empty;
            Q72Label.Content = string.Empty;
            M72Label.Content = string.Empty;
            Q73Label.Content = string.Empty;
            M73Label.Content = string.Empty;

            A81Label.Content = string.Empty;
            B81Label.Content = string.Empty;
            A82Label.Content = string.Empty;
            B82Label.Content = string.Empty;
            A83Label.Content = string.Empty;
            B83Label.Content = string.Empty;
            
            ALabel.Content = string.Empty;
            BLabel.Content = string.Empty;
        }

        // Analyse: normal or shorted normal
        private bool Analys()
        {
            //if count < 3
            int count = _resultCollection.Count;
            
            var mas = new double[count];
            for (int i = 0; i < count; i++)
            {
                mas[i] = _resultCollection[i].Ni;
            }
            
            var mas1 = new double[count / 3];
            for (int i = 0; i < count / 3; i++)
            {
                mas1[i] = mas[i];
            }

            var mas2 = new double[count / 3];
            for (int i = 0; i < count / 3; i++)
            {
                mas2[i] = mas[i + count / 3];
            }

            return mas1.Sum() > mas2.Sum();
        }
    }
}
