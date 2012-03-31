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
        Double dt;
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
            Double.TryParse(dtTextBox.Text, out dt);
            if (dt == 0)
                return;
            Int32.TryParse(KTextBox.Text, out K);
            if (K < 10)
                return;
            P = new Double[K];
            for(int i = 0; i < K; i++)
            {
                P[i] = CountProbability(i);
            }
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
            for (int i = 0; i < idx; i++)
            {
                Ni -= ResCol[i].ni;
            }
            return Ni;
        }

        double CountProbability(int i)
        {
            int i1 = i;
            int i2 = i + 1;
            return (double)(CountNi(i1) + CountNi(i2))/(2*N0);
        }

        void DrawExpGraphics()
        {
            ZedGraphControl ExpZedGraph = ExpWFH.Child as ZedGraphControl;
            GraphPane pane = ExpZedGraph.GraphPane;
            pane.CurveList.Clear();
            PointPairList list = new PointPairList();
            // Заполняем список точек
            for (int i = 0; i < K; i++)
            {
                list.Add( (ResCol[i].i - 0.5) * dt, P[i]);
            }
            LineItem myCurve = pane.AddCurve("Практически", list, System.Drawing.Color.Red, SymbolType.None);
            ExpZedGraph.AxisChange();
            // Обновляем график
            ExpZedGraph.Invalidate();
        }

        void DrawErlGraphics()
        {
            ZedGraphControl ErlZedGraph = ErlWFH.Child as ZedGraphControl;
            GraphPane pane = ErlZedGraph.GraphPane;
            pane.CurveList.Clear();
            PointPairList list = new PointPairList();
            // Заполняем список точек
            for (int i = 0; i < K; i++)
            {
                list.Add((ResCol[i].i - 0.5) * dt, P[i]);
            }
            LineItem myCurve = pane.AddCurve("Практически", list, System.Drawing.Color.Red, SymbolType.None);
            ErlZedGraph.AxisChange();
            // Обновляем график
            ErlZedGraph.Invalidate();
        }

        void DrawRelGraphics()
        {
            ZedGraphControl RelZedGraph = RelWFH.Child as ZedGraphControl;
            GraphPane pane = RelZedGraph.GraphPane;
            pane.CurveList.Clear();
            PointPairList list = new PointPairList();
            // Заполняем список точек
            for (int i = 0; i < K; i++)
            {
                list.Add((ResCol[i].i - 0.5) * dt, P[i]);
            }
            LineItem myCurve = pane.AddCurve("Практически", list, System.Drawing.Color.Red, SymbolType.None);
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
            LineItem myCurve = pane.AddCurve("Практически", list, System.Drawing.Color.Red, SymbolType.None);
            VejbZedGraph.AxisChange();
            // Обновляем график
            VejbZedGraph.Invalidate();
        }

        void DrawNormGraphics()
        {
            ZedGraphControl NormZedGraph = NormWFH.Child as ZedGraphControl;
            GraphPane pane = NormZedGraph.GraphPane;
            pane.CurveList.Clear();
            PointPairList list = new PointPairList();
            // Заполняем список точек
            for (int i = 0; i < K; i++)
            {
                list.Add((ResCol[i].i - 0.5) * dt, P[i]);
            }
            LineItem myCurve = pane.AddCurve("Практически", list, System.Drawing.Color.Red, SymbolType.None);
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
            LineItem myCurve = pane.AddCurve("Практически", list, System.Drawing.Color.Red, SymbolType.None);
            ShortNormZedGraph.AxisChange();
            // Обновляем график
            ShortNormZedGraph.Invalidate();
        }
    }
}
