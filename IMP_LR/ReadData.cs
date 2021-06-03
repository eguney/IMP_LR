using System;
using System.Collections.Generic;
using System.Text;

namespace IMP_LR
{
    public class ReadData
    {
        public dynamic Key { get; set; } = new { };
        public double W { get; set; }
        public int Count { get; set; }
        public double Prop { get; set; } = new Random().NextDouble();
        public double GetProbablity()
        {
            Prop = new Random().NextDouble();
            return Prop;
        }
    }
}
