using System;
using System.Collections.Generic;
using System.Text;

namespace IMP_LR
{
    public class BaseModel
    {
        public int Head { get; set; }
        public int Tail { get; set; }
    }

    public class BaseUInt16
    {
        public UInt16 Head { get; set; }
        public UInt16 Tail { get; set; }
    }

    public class ArcsUInt16Weighted
    {
        public UInt16 Head { get; set; }
        public UInt16 Tail { get; set; }
        public double Weight { get; set; }
    }
}
