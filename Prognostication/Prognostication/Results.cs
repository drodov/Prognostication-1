using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Prognostication
{
    class Results
    {
        public int i { set; get; }
        public int N { set; get; }
        public Results(int ii, int NN)
        {
            i = ii;
            N = NN;
        }
        public Results(int ii)
        {
            i = ii;
        }
    }
}
