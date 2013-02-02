namespace Prognostication
{
    class Results
    {
        public int I { set; get; }
        public double Ni { set; get; }
        
        public Results(int i, double ni)
        {
            I = i;
            Ni = ni;
        }

        public Results(int i)
        {
            I = i;
        }
    }
}
