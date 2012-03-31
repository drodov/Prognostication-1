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
            ErlZedGraph.GraphPane.Title.Text = "Закон Эрланга";
            RelZedGraph.GraphPane.Title.Text = "Закон Рэлея";
            VejbZedGraph.GraphPane.Title.Text = "Закон Вейбулла";
            NormZedGraph.GraphPane.Title.Text = "Нормальный закон";
            ShortNormZedGraph.GraphPane.Title.Text = "Усеченный нормальный закон";

            /*
            GraphPane pane = zedGraph.GraphPane;

            pane.CurveList.Clear();
            pane.Title.Text = "dfdf";

            PointPairList list = new PointPairList();

            double xmin = -50;
            double xmax = 50;

            // Заполняем список точек
            for (double x = xmin; x <= xmax; x += 0.01)
            {
                list.Add(x, f(x));
            }
            LineItem myCurve = pane.AddCurve("Sinc", list, System.Drawing.Color.Blue, SymbolType.None);

            zedGraph.AxisChange();

            // !!!
            // Линию рисуем после обновления осей с помощью AxisChange (), 
            // так как мы будем использовать значения
            // Нарисуем горизонтальную пунктирную линию от левого края до правого на уровне y = 0.5
            double level = 0.5;
            LineObj line = new LineObj(pane.XAxis.Scale.Min, level, pane.XAxis.Scale.Max, level);

            // Стиль линии - пунктирная
            line.Line.Style = System.Drawing.Drawing2D.DashStyle.Dash;

            // Добавим линию в список отображаемых объектов
            pane.GraphObjList.Add(line);

            // Нарисуем стрелку, указыающую на максимум функции
            // Координаты точки, куда указывает стрелка
            // Координаты привязаны к осям
            double xend = 0.0;
            double yend = f(0);

            // Координаты точки начала стрелки
            double xstart = xend + 5.0;
            double ystart = yend + 0.1;

            // Рисование стрелки с текстом
            // Создадим стрелку
            ArrowObj arrow = new ArrowObj(xstart, ystart, xend, yend);

            // Добавим стрелку в список отображаемых объектов
            pane.GraphObjList.Add(arrow);

            // Напишем текст около начала стрелки
            // Координаты привязаны к осям
            TextObj text = new TextObj("Max", xstart, ystart);

            // Отключим рамку вокруг текста
            text.FontSpec.Border.IsVisible = false;

            // Добавим текст в список отображаемых объектов
            pane.GraphObjList.Add(text);

            // Обновляем график
            zedGraph.Invalidate();*/

            dataGrid1.ItemsSource = ResCol;
        }

        private void KTextBox_TextChanged(object sender, TextChangedEventArgs e)
        {
            Int32 K;
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
    }
}
