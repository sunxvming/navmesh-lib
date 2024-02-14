using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace PathFinder
{
    class Program
    {
        const string PATHDLL = "navmesh";

        [DllImport(PATHDLL, CallingConvention = CallingConvention.Cdecl)]
        public static unsafe extern IntPtr CreatePath(double[] cont, int count);//int数组直接传
        [DllImport(PATHDLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern void SavePath(IntPtr handle, string resname);//字符串用byte[]
        [DllImport(PATHDLL, CallingConvention = CallingConvention.Cdecl)]
        public static unsafe extern IntPtr LoadPath(string resname);//字符串用byte[]
        [DllImport(PATHDLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern void FreePath(IntPtr handle);
        [DllImport(PATHDLL, CallingConvention = CallingConvention.Cdecl)]
        public unsafe static extern double* FindPath(IntPtr handle, double startx, double starty, double endx, double endy, out int size);
        [DllImport(PATHDLL, CallingConvention = CallingConvention.Cdecl)]
        public unsafe static extern double* FindCross(IntPtr handle, double startx, double starty, double facex, double facey);
        [DllImport(PATHDLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern int CheckPath(IntPtr handle, double startx, double starty, double endx, double endy);
        [DllImport(PATHDLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern IntPtr echoString(byte[] num);
        [DllImport(PATHDLL, CallingConvention = CallingConvention.Cdecl)]
        public unsafe extern static double* GetPoints(IntPtr cont, out int length);
        [DllImport(PATHDLL, CallingConvention = CallingConvention.Cdecl)]
        public unsafe static extern int* GetIndexs(IntPtr cont, out int length);


        static void Main(string[] args)
        {
            IntPtr ss = echoString(Encoding.ASCII.GetBytes("test string"));//传入返回字符串示例
            Console.WriteLine(Marshal.PtrToStringAnsi(ss));

            string path = Path.GetFullPath(".") + "/points.txt";
            unsafe
            {
                StreamReader reader = new StreamReader(path);
                List<double> numlist = new List<double>();
                string nextLine;
                while ((nextLine = reader.ReadLine()) != null)
                {
                    var a = Convert.ToDouble(nextLine);
                    numlist.Add(a);
                }
                reader.Close();

                double[] nums = numlist.ToArray();
        
                IntPtr navPath = CreatePath(nums, nums.Length);

                int length = 0;
                string pointpath = Path.GetFullPath(".") + "/points1.txt";
                StreamWriter stream = new StreamWriter(pointpath);
                double* cpoints = GetPoints(navPath, out length);
                for (int i = 0; i < length; i++)
                {
                    stream.WriteLine(cpoints[i]);
                }
                stream.Flush();
                stream.Close();


                int findSize = 0;
                double* finds = FindPath(navPath, 22f, 17.66f, 27.51f, 30.81f, out findSize);
                Console.WriteLine("----是否可行走--" + (findSize > 0 ? "是" : "否") + findSize);
                if (findSize > 0) {
                    for (int i = 0; i < findSize; i++)
                    {
                        Console.WriteLine(finds[i]);
                    }
                }

            }

        }

    }
}
