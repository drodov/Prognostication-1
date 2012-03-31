using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Prognostication
{
    class Results
    {
        public int i { set; get; }
        public int ni { set; get; }
        public Results(int ii, int Ni)
        {
            i = ii;
            ni = Ni;
        }
        public Results(int ii)
        {
            i = ii;
        }
    }
}
