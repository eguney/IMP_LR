using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;g
using System.Diagnostics;
//using ILOG.Concert;
//using ILOG.CPLEX;
//using System.Data.SqlClient;
//using Gurobi;
using System.Threading.Tasks;
using IMP_LR;
using System.Threading;

namespace IMP_LR
{
    class Program
    {
        static int N, N2, E, BIGM, trial_limit, k_init, seed, solution_method, SAA_M, SAA_N1, SAA_T, SAA_N2, SAA_N3, total_network_depth, Cplex_seconds, k_max, run_pipage, run_CELF, run_LR, run_MIP, run_LATIN, run_reduction, run_MEM, theoption, runLT, runTIPTOP, LR_NG_Cuts, no_of_LR_iter;
        static int symmetric_network, save_model, save_log2, mip_model_version, max_network_depth, mip_with_cuts, networkID, is_BUDGETED, counter_prob, counter_prob2;
        static uint m_w, m_z;
        static int[] no_of_active_set, indegree, outdegree, SAA_total_network_depth, SAA_max_network_depth, memory_sample;
        static long[] node_score;
        static double prop_prob, SA_1_hat_obj, SA_3_obj, param_epgap, R1_avg, R1_var, R1_std, T_avg, T_var, T_std, T2_avg, T2_var, T2_std, prob_scaler, pipage_objective_LP, pipage_objective_IP, pipage_T, tiptop_eps, tiptop_rho;
        static double CELF_z, CELF_t, z_pipageInf, z_LRUB, z_LRLB, zLRinf, t_LR, fixed_probability, mem_limit;
        static double[] CELF_z_arr, CELF_t_arr, y_obj;
        static int[,] arcs_int;
        static int[] node_index;
        static int[,] arcs_intID;
        static bool[,] x_exist;
        static int[,,] network_depth;
        static string filename, pipage_solution, CELF_solution, LR_1_solution, pipage_count, solution_LR;
        static string[] file_names, CELF_sol_arr;
        static string[,,] SAA_neigh, SAA_pred, SAA_tree, SAA_fwd_tree;
        static double[] avg_no_AS, node_threshold, SA_1_obj, SA_1_test_obj, SA_2_obj, recruitment_cost, T_LB, T2_LB, y_i, LP_duals;
        static double[,] influence_prob, SAA_1_prb, SAA_2_prb, SAA_3_prob, Beta_ir, x_ir;
        static float[,] SAA_float;
        static double[,,] SAA_prob, SAA_2_prob;
        static List<int> my_node_set, node_set, init_set, CA_set, active_set, NA_set, selected_set, result_set, SAA_current_list, selected_set2;
        static List<double> arc_weight, arc_distance, accumulated_influence, weight_set;
        static List<string> neigh_set, arc_set, arcID_set, pred_set, neigh_edge_set;
        static List<List<int>> neigh_list, pred_list, Tiptop_pred_list, hypernodes;
        static List<List<int>> neigh_edge_list;
        static List<List<int>> SAA_list, SAA_pipage_list, SAA_LR_list;
        static List<List<double>> SAA_3_prob_dbl, Tiptop_w_list;
        static List<List<List<int>>> SAA_neigh_list, SAA_pred_list, SAA_tree_list, SAA_fwd_tree_list;
        static StreamWriter swvar, swvar_summary;
        //static Cplex cplex_model;
        //static ILPMatrix lp_model;
        static Stopwatch sw;
        static TimeSpan sw_elapsed;
        static StreamWriter swtemp2;
        static List<ReadData> unique_arcs;
        static List<ArcsUInt16Weighted> only_arcs;

        static void Main(string[] args)
        {
            //swvar_summary=new StreamWriter("results.txt");
            sw = new Stopwatch();

            file_names = new string[10] { "arcs2000", "phy-500", "phy-1k", "phy-2k", "phy-5k", "phy-10k", "phy-20k", "phy-50k", "phy-100k", "phy-200k" };
            //file_names = new string[7] { "run_1K", "run_2K", "run_5K", "run_10K", "run_20K", "run_50K", "run_100K" };
                string[] file_names3 = new string[6] { "Email-Enron.im", "phy.im", "p2p-Gnutella04.im", "CollegeMsg_comma.txt", "Slashdot0902.txt", "test.txt" };

            int[] SAA_N1_arr = new int[16] { 5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200 };
            int[] fixed_k = new int[7] { 2, 3, 4, 5, 10, 15, 25 };
            int[] fixed_N1 = new int[4] { 30, 200, 300, 500 };
            double[] fixed_probs = new double[2] { 0.01, 0.05 };
            //swvar = new StreamWriter("log-degreecentrality.txt");
            SAA_M = 1;
            SAA_N1 = 200;
            SAA_N2 = 5000;
            SAA_N3 = 5000;
            SAA_T = 1;
            prob_scaler = 1;
            run_LR = 1;
            LR_NG_Cuts = 0;
            run_LATIN = 0;
            run_reduction = 0;
            runLT = 0;
            int file_index_last = file_names.Length;
            file_index_last = 5;
            k_max = 7;

            int file_index = 3;
            filename = file_names3[file_index];
            //filename = file_names[file_index];

            fixed_probability = 0.25;


            initialize_sinan();
           
            
            // Main for loop for different k values (seed set sizes)
            for (int k = 4; k < k_max; k = k + 1) //kmax=5
            {
                k_init = fixed_k[k];
                swvar = new StreamWriter("log-" + filename + "_" + k + ".txt");
                
                //BUDGET = k_init;
                trial_limit = SAA_N1;  // number of scenarios/samples

                System.Console.WriteLine("Finished Initialization");

                System.Console.WriteLine("Iteration K:" + k + "  data: " + filename + " R:" + SAA_N1 + " BUDGET:" + k_init + " p:" + fixed_probability);


                sw.Start(); // Stopwatch start

                SAA_Main_Loop(); // data is read! Now prepare the RR sets and then solve it.

                sw.Reset();

                swvar.Close();
            }


            System.Console.WriteLine("File " + filename + " finished");
        }


        public static void initialize_sinan()   // read data, create nodes/arcs/B/c_i
        {
            symmetric_network = 0; // 0-> directed network, 1-> undirected network
            int arc_weights_ON = 0;
            int header_row = 0;
            int fix_prob = 0;

            m_w = 521288629;
            m_z = 362436069;
            BIGM = 1000000;
            //SetSeedFromSystemTime();

            // read raw data from file
            // Please correct this location for your text file

            string path;
            //path = @"C:\Users\Evren\Downloads\inst\" + filename;
            //path = @"C:\data\CollegeMsg.im";
            path = @"C:\data\"+filename;

            arc_set = new List<string>();

            var readData = (File.ReadLines(path).Skip(0).Select(a => {
                return new { Id = a };
            }).AsParallel().ToList());

            if (arc_weights_ON == 0)
            {
                unique_arcs = readData.GroupBy(a => a.Id).Select(a => new { Key = a, Count = a.Count() }).Select(a =>
            new ReadData
            {
                Key = new { Head = int.Parse(a.Key.Key.Split(',')[0]), Tail = int.Parse(a.Key.Key.Split(',')[1]) },
                Count = a.Count,
                W = (a.Count == 1) ? fixed_probability : (1 - Math.Pow((1 - fixed_probability), a.Count)),
            }).AsParallel().ToList();
            }

            else
            {
                unique_arcs = readData.GroupBy(a => a.Id).Select(a => new { Key = a, Count = a.Count() }).Select(a =>
                new ReadData
                {
                    Key = new { Head = int.Parse(a.Key.Key.Split(',')[0]), Tail = int.Parse(a.Key.Key.Split(',')[1]) },
                    Count = a.Count,
                    W = double.Parse(a.Key.Key.Split(',')[2]),
                }).AsParallel().ToList();
            }


            //new ReadData
            //{
            //    Key = new { Head = int.Parse(a.Key.Key.Split(',')[0]), Tail = int.Parse(a.Key.Key.Split(',')[1]) },
            //    Count = a.Count,
            //    W = Double.Parse(a.Key.Key.Split(',')[2]),
            //}).AsParallel().ToList();

            var nodes = unique_arcs.Select(a => a.Key.Head).Distinct().Union(unique_arcs.Select(b => b.Key.Tail).Distinct()).Distinct().AsParallel().ToList();
            //var left_nodes_neighbours = unique_arcs.GroupBy(a => a.Key.Head).Select(a => new { Head = a.Key, RefList = a.Select(b => b.Key.Tail).ToList() }).AsParallel().ToList();


            //List<ReadData> UniqueArcs = new List<ReadData>();
            //UniqueArcs.AddRange(unique_arcs);
            //UniqueArcs = unique_arcs.Select(a => a).AsParallel().ToList();



            result_set = new List<int>();
            N = nodes.Count;
            E = unique_arcs.Count;
            node_set = new List<int>();
            for (int i = 0; i < N; i++)
                node_set.Add((int)nodes[i]);


            System.Console.WriteLine("Finished Advanced Initialization");

        }

        public static void initialize_SAA_neighbourhood_fast(int sample_size, List<dynamic> reslist)
        {
            System.Console.WriteLine("Starting Fast SAA_Neigh");
            SAA_neigh_list = new List<List<List<int>>>();
            SAA_pred_list = new List<List<List<int>>>();
            int counter2 = 0;
            x_exist = new bool[sample_size, N];
            //List<DirectedListGraph<int,UnweightedEdge>> m_list = new List<DirectedListGraph<int, UnweightedEdge>>();

            int Nmax = node_set.Max() + 1;
            node_index = new int[Nmax];
            for (int i = 0; i < N; i++)
            {
                node_index[node_set[i]] = i;
            }
            //x_exist = new bool[N, sample_size];

            for (int r = 0; r < sample_size; r++)
            {
                SAA_neigh_list.Add(new List<List<int>>());
                SAA_pred_list.Add(new List<List<int>>());
                //var m = new DirectedListGraph<int, UnweightedEdge>();
                for (int i = 0; i < N; i++)
                {
                    //x_exist[i, r] = true;
                    SAA_neigh_list[r].Add(new List<int>());
                    SAA_pred_list[r].Add(new List<int>());
                    //m.AddVertex(i);
                }
                //m_list.Add(m);
            }

            int head, tail;

            for (int r = 0; r < sample_size; r++)
            //Parallel.For(0, sample_size, r =>
            {
                foreach (BaseModel item in reslist[r])
                {
                    head = node_index[(int)item.Head];
                    tail = node_index[(int)item.Tail];
                    SAA_neigh_list[r][head].Add(tail);
                    SAA_pred_list[r][tail].Add(head);
                    counter2++;
                    //m_list[r].AddEdge(new UnweightedEdge(head, tail));
                }
            }
            //);
            //swtemp.Close();
            System.Console.WriteLine("Finished Fast Neigh! ...... Counter:" + E * SAA_N1 + " & Counter2:" + counter2);
            //determine_network_trees_IM(sample_size, m_list);
        }  //create stochastic neighbourhood list for each node, for each sample, for each batch


        public static void initialize_neighbourhood_fast_Tiptop()
        {
            System.Console.WriteLine("Starting Edge Neighbourhood for Tiptop");
            Tiptop_w_list = new List<List<Double>>();
            Tiptop_pred_list = new List<List<int>>();
            int counter = 0;

            int Nmax = node_set.Max() + 1;
            node_index = new int[Nmax];
            for (int i = 0; i < N; i++)
            {
                node_index[node_set[i]] = i;
                Tiptop_pred_list.Add(new List<int>());
                Tiptop_w_list.Add(new List<double>());
            }
            ////x_exist = new bool[N, sample_size];

            int head, tail;
            Double weight;

            foreach (var item in unique_arcs)
            {
                head = node_index[(int)item.Key.Head];
                tail = node_index[(int)item.Key.Tail];
                weight = item.W;
                Tiptop_w_list[tail].Add(weight);
                Tiptop_pred_list[tail].Add(head);
                counter++;
            }
            //}
            ////);
            ////swtemp.Close();
            System.Console.WriteLine("Finished Fast Neigh! ...... Counter:" + counter);
            ////determine_network_trees_IM(sample_size, m_list);
        }


        public static void initialize_random_arcs()
        {
            var reslist = new List<dynamic>();
            var tasks = new List<Task<List<BaseModel>>>();
            var rand = new Random();
            for (int i = 0; i < E; i++)
            {
                unique_arcs[i].Prop = rand.NextDouble();
            }

            for (int i = 0; i < SAA_N1; i++)
            {
                var luckyarcs = Task.Run<List<BaseModel>>(() => unique_arcs.Select(a => a).Where(a => (a.W >= a.GetProbablity())).AsParallel().ToList().Select(b =>
                            new BaseModel { Head = b.Key.Head, Tail = b.Key.Tail }).AsParallel().ToList());
                tasks.Add(luckyarcs);
                //var luckyarcs = unique_arcs.Select(a => a).Where(a => (a.W >= GetUniform())).AsParallel().ToList().Select(b =>
                //               new BaseModel { Head = b.Key.Head, Tail = b.Key.Tail }).AsParallel().ToList();
                //reslist.Add(luckyarcs);
            }
            var results = Task.WhenAll(tasks);
            foreach (var task in tasks)
            {
                var l1 = task.Result;
                reslist.Add(l1);
            }
            initialize_SAA_neighbourhood_fast(SAA_N1, reslist);
        }

        public static void initialize_random_arcs_LT(int sample_size)
        {
            System.Console.WriteLine("Starting Fast SAA_Neigh LT");
            SAA_neigh_list = new List<List<List<int>>>();
            SAA_pred_list = new List<List<List<int>>>();
            int counter2 = 0;
            x_exist = new bool[sample_size, N];

            int Nmax = node_set.Max() + 1;
            node_index = new int[Nmax];
            for (int i = 0; i < N; i++)
            {
                node_index[node_set[i]] = i;
            }
            //x_exist = new bool[N, sample_size];

            for (int r = 0; r < sample_size; r++)
            {
                SAA_neigh_list.Add(new List<List<int>>());
                SAA_pred_list.Add(new List<List<int>>());
                //var m = new DirectedListGraph<int, UnweightedEdge>();
                for (int i = 0; i < N; i++)
                {
                    //x_exist[i, r] = true;
                    SAA_neigh_list[r].Add(new List<int>());
                    SAA_pred_list[r].Add(new List<int>());
                    //m.AddVertex(i);
                }
                //m_list.Add(m);
            }

            neigh_list = new List<List<int>>();
            for (int i = 0; i < N; i++)
            {
                neigh_list.Add(new List<int>());
            }

            int current_head = 0; ;
            int current_tail = 0; ;

            for (int i = 0; i < E; i++)
            {
                current_head = node_index[(int)unique_arcs[i].Key.Head];
                current_tail = node_index[(int)unique_arcs[i].Key.Tail];
                neigh_list[current_head].Add(current_tail);
            }

            var rand = new Random();
            double rand_number = 0;
            int neigh_order = 0;

            for (int r = 0; r < sample_size; r++)
            {
                for (int i = 0; i < N; i++)
                {
                    neigh_order = neigh_list[i].Count;
                    if (neigh_order > 0)
                    {
                        rand_number = Math.Floor(rand.NextDouble() * neigh_order);
                        SAA_neigh_list[r][i].Add(neigh_list[i][(int)rand_number]);
                        SAA_pred_list[r][(int)neigh_list[i][(int)rand_number]].Add((int)i);
                    }
                }
            }
        }

        public static List<double> Shuffle(List<double> array)
        {
            int n = array.Count;
            Random rng = new Random();
            double tmp = 0;
            while (n > 1)
            {
                int k = rng.Next(n--);
                tmp = array[n];
                array[n] = array[k];
                array[k] = tmp;
            }
            return array;
        }

        public static void initialize_SAA_Float()  // create probabilities for objective function calculation diffusion with a large sample size
        {
            float prob;
            SAA_float = new float[E, SAA_N2];

            for (int c = 0; c < E; c++)
            {
                for (int b = 0; b < SAA_N2; b++)
                {
                    //for (int a = 0; a < SAA_M; a++)
                    {
                        prob = (float)GetUniform();
                        SAA_float[c, b] = prob;
                    }
                }
            }
        }
        
        public static void SAA_Main_Loop()
        {
            y_obj = new double[N];  // keeps the increments after removing singletons

            double max_obj = -1;
            
            SAA_LR_list = new List<List<int>>(SAA_M);
            string SAA_str_solution = "";
            Double dt_LR = 0;

            double[] LR_UB_arr = new double[SAA_M];
            double[] LR_LB_arr = new double[SAA_M];
            double[] LR_inf_arr = new double[SAA_M];
            string[] LR_sol_arr = new string[SAA_M];
            
            double LR_tot_UB = 0; double LR_tot_LB = 0; double LR_avg_UB = 0; double LR_avg_LB = 0; double LR_best_inf = 0; int LR_best_index = -1; string LR_best_seed = "";
            double LR_UB_std = 0;
            double LR_inf_std = 0;
            double LR_total_inf = 0;


            if ( run_LR == 1)
            {
                Stopwatch sw_LR = new Stopwatch();
                Stopwatch swx = new Stopwatch();
                double temp_obj = 0;
                int max_index = -1;
                swx.Start();
                for (int i = 0; i < SAA_M; i++)
                {

                    if (runLT == 1)
                        initialize_random_arcs_LT(SAA_N1);
                    else
                        initialize_random_arcs();  // live-arc operation

                    determine_network_trees(SAA_N1); // construct RR sets

                    swx.Stop(); System.Console.WriteLine(swx.Elapsed.Ticks); swx.Reset();

                    if (run_reduction == 1)   // if you want to remove singletons
                        reduce_tree_list();

                    sw_LR.Start();
                    //z_LRUB is written in function, z_LRLB is what function returns (feasible solution), zLRinf, t_LR
                    // Here is the Lagrangean relaxation with subgradient optimization subroutine!!!

                    LR_LB_arr[i] = (double)LR_0(i, 0, SAA_N1) / SAA_N1;  //feasible solution's objective value is a lowerbound

                    LR_UB_arr[i] = z_LRUB;  // LR functions final value is an upper bound (since it is a relaxation)
                    LR_sol_arr[i] = solution_LR;  // seed set
                    sw_LR.Stop();
                    double sw_LR_elapsed = sw_LR.ElapsedMilliseconds; t_LR = sw_LR_elapsed;

                }


                //reporting
                //LR results
                // The seed set is put in simulation for an actual influence estimation. Sample size is much much larger

                if (run_LR == 1)
                {
                    int inf_method = 0;
                    if (SAA_M > 100)
                    {
                        initialize_SAA_Float();
                        inf_method = 6;
                    }
                    else
                    {
                        inf_method = 2;
                    }
                    for (int i = 0; i < SAA_M; i++)
                    {

                        LR_inf_arr[i] = evaluate_influence(SAA_LR_list[i], i, SAA_N2, inf_method);
                        //LR_inf_arr[i] = i+1;
                        System.Console.WriteLine(i + "-->" + LR_UB_arr[i] + "/" + LR_LB_arr[i] + "  inf: " + LR_inf_arr[i]);
                        if (LR_inf_arr[i] > LR_best_inf)
                        {
                            LR_best_seed = LR_sol_arr[i]; LR_best_inf = LR_inf_arr[i]; LR_best_index = i; System.Console.WriteLine("New best found at batch " + (i + 1));
                        }
                        LR_tot_UB = LR_tot_UB + LR_UB_arr[i]; LR_tot_LB = LR_tot_LB + LR_LB_arr[i]; LR_total_inf = LR_total_inf + LR_inf_arr[i];
                    }
                    LR_avg_UB = LR_tot_UB / SAA_M; LR_avg_LB = LR_tot_LB / SAA_M;

                    for (int i = 0; i < SAA_M; i++)
                    {
                        LR_UB_std = LR_UB_std + (LR_UB_arr[i] - LR_avg_UB) * (LR_UB_arr[i] - LR_tot_UB / SAA_M);
                        LR_inf_std = LR_inf_std + (LR_inf_arr[i] - LR_total_inf / SAA_M) * (LR_inf_arr[i] - LR_total_inf / SAA_M);
                    }
                    LR_UB_std = Math.Sqrt(LR_UB_std / SAA_M); LR_inf_std = Math.Sqrt(LR_inf_std / SAA_M);


                    //sw.Stop(); sw_elapsed = sw.Elapsed; pipage_T = System.Convert.ToDouble(sw_elapsed.Ticks);
                    System.Console.WriteLine("\r\n ---- LR-1 Solutions---");
                    System.Console.WriteLine("Sample Averages (IP/LP): " + LR_avg_UB + "/" + LR_avg_LB);
                    System.Console.WriteLine("Best seed: " + LR_sol_arr[LR_best_index]);
                    System.Console.WriteLine("Best Objective Value Estimate: " + LR_inf_arr[LR_best_index]);
                    System.Console.WriteLine("Stdevs: " + LR_UB_std + " - " + LR_inf_std);
                    System.Console.WriteLine("Time : " + System.Convert.ToDouble(t_LR / 10000000.0));
                }

                sw.Reset(); sw_LR.Reset();
            }


            string sqltext;
            //SqlConnection conn = new SqlConnection(constr);
            //conn.Open();
            //SqlCommand command = conn.CreateCommand();

            sqltext = " VALUES (" + SAA_M + "," + SAA_N1 + "," + SAA_N2 + "," + SAA_N3 + ",1,1," + k_init + ",1," + N + "," + E + ",-1," + SA_1_hat_obj + "," + max_obj + "," + SA_3_obj + ",'" + SAA_str_solution + "'," + prop_prob + ",'" + k_init + "','" + R1_std + "','" + T2_std + "'," + pipage_objective_LP + "," + pipage_objective_IP + "," + z_pipageInf + "," + pipage_T + ",'" + pipage_solution + "','" + pipage_count + "'," + LR_avg_UB + "," + LR_avg_LB + "," + LR_best_inf + "," + t_LR + ",'" + LR_best_seed + "'," + LR_UB_std + "," + LR_inf_std + "," + fixed_probability + "," + no_of_LR_iter + ")";
            sqltext = sqltext.Replace("'NaN'", "'0'");


            //command.Dispose();
            //conn.Close();
            SAA_tree_list.Clear();
            //SAA_fwd_tree_list.Clear();
        }

        public static string fwd_tree_with_BFS(int starting_node, int r, int m) // determine the reachable set by Breadth-first-search for each node - sample couple
        {
            List<int> temp_nodes = new List<int>();
            List<int> visited_nodes = new List<int>();
            int[] node_colors = new int[N];
            int[] node_labels = new int[N];
            int starting_node_ID = node_set[starting_node];
            for (int i = 0; i < N; i++)
            {
                node_colors[i] = 0;
                node_labels[i] = -1;
            }
            int node_index = -1;
            string[] neigh_arr;
            int curr_index = 0;

            temp_nodes.Add(starting_node);
            //visited_nodes.Add(starting_node);

            node_index = starting_node;
            int statint_node_index = starting_node;
            node_colors[node_index] = 1;
            node_labels[node_index] = 0;

            while (temp_nodes.Count != 0)
            {
                node_index = temp_nodes[0];
                if (SAA_neigh[node_index, r, m] != null)
                {
                    neigh_arr = SAA_neigh[node_index, r, m].Split(new char[] { ',' });     //neighbour stringini split edelim
                    foreach (string str in neigh_arr)   //tüm komşuları aktifleştirmeye çalış
                    {
                        curr_index = Convert.ToUInt16(str);
                        if (node_colors[curr_index] == 0)
                        {
                            temp_nodes.Add(curr_index);
                            if (!visited_nodes.Contains(curr_index))
                                visited_nodes.Add(curr_index);
                            node_colors[curr_index] = 1;
                            if (curr_index != statint_node_index && node_labels[curr_index] == -1)
                                node_labels[curr_index] = node_labels[node_index] + 1;
                        }

                    }
                }
                temp_nodes.Remove(temp_nodes[0]);
                node_colors[node_index] = 2;
            }
            StringBuilder builder = new StringBuilder();
            foreach (int node in visited_nodes)
            {
                // Append each int to the StringBuilder overload.
                builder.Append(node).Append(",");
            }
            string result = builder.ToString();
            result = result.TrimEnd(',');

            return result;
        }

        //new version
        public static void determine_network_trees(int sample_size) // Determine the RR set by BFS main function
        {
            System.Console.WriteLine("Reverse Reachable set for : " + SAA_M * sample_size * N + " node-run-sample couples (takes 1-2 minutes)");
            SAA_tree_list = new List<List<List<int>>>();

            int[] tree_size = new int[N];
            node_score = new long[N];
            int[] tree = new int[N];
            //for (int m = 0; m < SAA_M; m++)
            {
                int m = 0;
                //  System.Console.WriteLine("m:" + m);
                for (int r = 0; r < sample_size; r++)
                {
                    if (r % 5 == 0)
                    {// System.Console.Write("-" + r); System.Console.WriteLine(" | Mem: " + mem + "MB , ts: "+ tree_size.Sum());
                        System.Console.WriteLine("---" + r);
                    }

                    {
                        SAA_tree_list.Add(new List<List<int>>());

                        for (int i = 0; i < N; i++)
                        {
                            SAA_tree_list[r].Add(new List<int>());
                            SAA_tree_list[r][i].Add(i);
                            node_score[i]++;

                            if (i % 10000 == 0)
                                System.Console.Write("." + i);
                            //add visited_nodes

                            if (SAA_pred_list[r][i].Count() > 0)
                            {
                                SAA_tree_list[r][i].AddRange(BFS_IM(i, r, m, tree, SAA_pred_list, SAA_tree_list));
                                tree[i] = 1;  // if node i has predecessors let's flag it to speed up BFS
                            }
                            else
                            {
                                y_obj[i]++;  // if this node does not have predecessor it is a singleton
                            }
                        }
                        for (int i = 0; i < N; i++)
                            tree[i] = 0;
                    }
                }

                //System.Console.WriteLine(m);
            }
            //Int64 totalsize = tree_size.Sum();
            System.Console.WriteLine("RR Set finished! Size: ");
            //determine_network_trees2(sample_size);

        }

        public static List<int> BFS_IM(int starting_node, int r, int m, int[] tree, List<List<List<int>>> neigh, List<List<List<int>>> treelist) // determine the reachable set by Breadth-first-search for each node - sample couple
        {
            var vertices = new Queue<int>();
            vertices.Enqueue(starting_node);

            HashSet<int> visited_nodes = new HashSet<int>();
            bool[] marked = new bool[N];
            while (vertices.Any())
            {
                int curVertex = vertices.Dequeue();
                marked[curVertex] = true;
                foreach (int vertex in neigh[r][curVertex].Where(x => !marked[x]))
                {
                    if (tree[vertex] == 1)   // If we already have the BFS list for this vertex directly add those
                    {
                        //node_labels[node] = 0;
                        foreach (int nodes in treelist[r][vertex])
                        {
                            visited_nodes.Add(nodes);
                            marked[nodes] = true;
                        }
                    }
                    else
                    {
                        visited_nodes.Add(vertex);
                        vertices.Enqueue(vertex);
                    }

                }
            }

            if (visited_nodes.Contains(starting_node))
            {
                visited_nodes.Remove(starting_node);
            }
            foreach (int node in visited_nodes)
                node_score[node]++;
            return visited_nodes.ToList();
        }

        public static List<int> BFS_IM_simple(int starting_node, int r, int[,] tree, List<List<List<int>>> neigh, List<List<List<int>>> treelist) // determine the reachable set by Breadth-first-search for each node - sample couple
        {
            var vertices = new Queue<int>();
            vertices.Enqueue(starting_node);

            HashSet<int> visited_nodes = new HashSet<int>();
            bool[] marked = new bool[N];
            while (vertices.Any())
            {
                int curVertex = vertices.Dequeue();
                marked[curVertex] = true;
                foreach (int vertex in neigh[r][curVertex].Where(x => !marked[x]))
                {
                    if (tree[r, vertex] == 1)
                    {
                        //node_labels[node] = 0;
                        foreach (int nodes in treelist[r][vertex])
                        {
                            visited_nodes.Add(nodes);
                            marked[nodes] = true;
                        }
                    }
                    else
                    {
                        visited_nodes.Add(vertex);
                        vertices.Enqueue(vertex);
                    }

                }
            }

            if (visited_nodes.Contains(starting_node))
            {
                visited_nodes.Remove(starting_node);
            }
            return visited_nodes.ToList();
        }

        public static List<int> BFS_IM_simple2(int starting_node, int r, int[] tree, List<List<List<int>>> neigh, List<List<int>> treelist) // determine the reachable set by Breadth-first-search for each node - sample couple
        {
            var vertices = new Queue<int>();
            vertices.Enqueue(starting_node);

            HashSet<int> visited_nodes = new HashSet<int>();
            bool[] marked = new bool[N];
            while (vertices.Any())
            {
                int curVertex = vertices.Dequeue();
                marked[curVertex] = true;
                foreach (int vertex in neigh[r][curVertex].Where(x => !marked[x]))
                {
                    if (tree[vertex] == 1)
                    {
                        //node_labels[node] = 0;
                        foreach (int nodes in treelist[vertex])
                        {
                            visited_nodes.Add(nodes);
                            marked[nodes] = true;
                        }
                    }
                    else
                    {
                        visited_nodes.Add(vertex);
                        vertices.Enqueue(vertex);
                    }

                }
            }

            if (visited_nodes.Contains(starting_node))
            {
                visited_nodes.Remove(starting_node);
            }
            return visited_nodes.ToList();
        }

        public static List<int> BFS2(int starting_node, int r, int m, List<List<List<int>>> neigh) // determine the reachable set by Breadth-first-search for each node - sample couple, this is  a simple on with no check on intermediate nodes RR-list
        {
            var vertices = new Queue<int>();
            vertices.Enqueue(starting_node);

            HashSet<int> visited_nodes = new HashSet<int>();
            bool[] marked = new bool[N];
            while (vertices.Any())
            {
                int curVertex = vertices.Dequeue();
                marked[curVertex] = true;
                foreach (int vertex in neigh[r][curVertex].Where(x => !marked[x]))
                {
                    visited_nodes.Add(vertex);
                    vertices.Enqueue(vertex);
                }
            }

            if (visited_nodes.Contains(starting_node))
            {
                visited_nodes.Remove(starting_node);
            }

            return visited_nodes.ToList();
        }

        public static List<int> BFS_TipTop(int starting_node) // BFS for Tiptop
        {
            var vertices = new Queue<int>();
            vertices.Enqueue(starting_node);
            int n_visit_mark = 0;
            Random rnd = new Random();

            HashSet<int> visited_nodes = new HashSet<int>();
            bool[] marked = new bool[N];
            int[] visit_mark = new int[N];
            visit_mark[n_visit_mark++] = starting_node;
            marked[starting_node] = true;
            //marked[curVertex] = true;
            while (vertices.Any())
            {
                int curVertex = vertices.Dequeue();

                for (int i = 0; i < Tiptop_pred_list[curVertex].Count; i++)
                {
                    int item = Tiptop_pred_list[curVertex][i];
                    if (marked[item])
                        continue;
                    else
                    {
                        double temp_arc_weight = rnd.NextDouble();
                        if (temp_arc_weight > Tiptop_w_list[curVertex][i])
                            continue;
                        visit_mark[n_visit_mark++] = item;
                        marked[item] = true;
                        visited_nodes.Add(item);
                    }
                    vertices.Enqueue(item);
                }
            }

            if (visited_nodes.Contains(starting_node))
            {
                visited_nodes.Remove(starting_node);
            }
            foreach (int node in visited_nodes)
                node_score[node]++;
            for (int i = 0; i < n_visit_mark; i++)
                marked[visit_mark[i]] = false;
            return visited_nodes.ToList();
        }

        public static void reduce_tree_list()
        {
            //determine_fwd_network_trees_list(SAA_N1);
            List<int> dominated = new List<int>();
            int count = 0;
            count = countlist(SAA_tree_list);
            int no_of_x_removed = 0;
            for (int i = 0; i < N; i++)
                y_obj[i] = 0;

            for (int i = 0; i < 1; i++)
            {
                //if (SAA_fwd_tree_list[0][i].Count > 10)
                goto gotoloc2;
                y_obj[i] = 0;
                List<int> PL2 = new List<int>();
                PL2.AddRange(SAA_tree_list[0][i]);  //take all predecessors of node i from the first sample
                PL2.Remove(i);
                foreach (int node in PL2.ToList())    // check each predecessor of i one by one, if it exists in pred list of i in all samples
                {
                    for (int r = 1; r < SAA_N1; r++)
                    {
                        //List<int> list_to_check = new List<int>();
                        //list_to_check.AddRange(SAA_fwd_tree_list[r][i]);
                        foreach (int location in SAA_fwd_tree_list[r][i].ToList())
                        {
                            if (SAA_tree_list[r][(int)location].IndexOf(node) < 0)
                            {
                                PL2.Remove(node);
                                //break;
                                goto gotoloc;
                            }
                        }
                    }
                    dominated.Add(i);
                    remove_i_from_list(i);
                    break;
                gotoloc:;
                }
            gotoloc2:;
            }
            System.Console.WriteLine("Removed: ");
            count = countlist(SAA_tree_list);
            int temp_removed = 0;

            for (int r = 0; r < trial_limit; r++)
            {
                for (int i = 0; i < N; i++)
                {
                    if (SAA_tree_list[r][i].Count == 1)
                    {
                        int thenode = (int)SAA_tree_list[r][i][0];
                        y_obj[thenode] = y_obj[thenode] + 1;
                        remove_constraint(i, thenode, r);
                        no_of_x_removed++;
                    }
                }           //end of for loop for nodes

                //System.Console.Write(no_of_x_removed-temp_removed+"-");
                temp_removed = no_of_x_removed;
            }
            System.Console.WriteLine("Dominated: " + dominated.Count);
            System.Console.WriteLine("no_of_x_removed: " + no_of_x_removed);

            //return dominated;
        }

        public static void remove_i_from_list(int node)
        {
            int thenode = (int)node;
            {
                for (int r = 0; r < SAA_N1; r++)
                {
                    List<int> list_to_check = new List<int>();
                    list_to_check.AddRange(SAA_fwd_tree_list[r][thenode]); // find the locations where node exists
                    foreach (int location in list_to_check.ToList())
                    {
                        SAA_tree_list[r][(int)location].Remove(thenode);    //RR setten sil
                    }
                    SAA_fwd_tree_list[r][thenode].Clear();  //Fwd setten sil
                }
            }
        }

        public static void remove_constraint(int i, int node, int r)
        {
            //double mem = 0;
            //int x=SAA_tree_list[r][i].Count;
            //mem = System.GC.GetTotalMemory(true) ;
            SAA_tree_list[r][i].Clear();
            //mem = System.GC.GetTotalMemory(true);
            //x = SAA_tree_list[r][i].Count;
            SAA_tree_list[r][i] = null;
            //mem = System.GC.GetTotalMemory(true);


            // SAA_fwd_tree_list[r][node].Remove(i);
            x_exist[i, r] = false;

        }

        public static int countlist(List<List<List<int>>> list)
        {
            int count = 0;
            for (int i = 0; i < N; i++)
            {
                for (int r = 0; r < SAA_N1; r++)
                {
                    count = count + list[r][i].Count;
                }
            }
            return count;
        }


        public static double LR_0(int SAA_m, int remove_nodes, int sample_size)
        {
            double solution = 0.0;
            no_of_LR_iter = 0;
            Beta_ir = new double[N, sample_size];   // Lagrange multipliers
            x_ir = new double[N, sample_size];
            y_i = new double[N];
            double[] c_i = new double[N];   // Coefficients for y-variables

            double z_UB = 50.0 * BIGM; //min UB
            double z_LB = -50.0 * BIGM;  // feas sol.
            double lambda = 2;              // Parameter used in subgradient optimization
            int[,] gradient = new int[N, sample_size];
            double temp_UB = 0;
            double temp_LB = 0;
            double step_size = 0;
            double stopping_threshold = 0.001;    // when UB and LB converge this much we stop
            int max_no_improvement_limit = 30;  // # ofsubgradient iterations before halving Lambda
            int LR_continue = 1;
            double z_UB_prev = 0;
            int no_improvement = 0;
            int iterlim = 1500;             // max # of subgradient iterations allowed
            int iter = 1;
            int[,] tree = new int[sample_size, N];
            List<int> y_list = new List<int>();
            List<int> y_iter = new List<int>();
            List<int> y_best = new List<int>();
            List<int> y_toremove = new List<int>();
            List<int> temp_y = new List<int>();
            List<int> y_ever = new List<int>();
            List<List<List<int>>> y_influencees = new List<List<List<int>>>();
            double t1, t2, t3, t4, t5, t6;
            sw.Reset();
        
            //initialize beta_ir
            for (int i = 0; i < N; i++)
            {
                y_iter.Add(i);

                c_i[i] = node_score[i];// + y_obj[i]; // Node scores are the initial C_i values calculated during BFS runs

                for (int r = 0; r < sample_size; r++)
                {
                    Beta_ir[i, r] = 1;  // we initialize all Lagrange multipliers to 1
                }
            }
            System.Console.WriteLine("Finished initialization of LR! SAA" + SAA_m);

            //sw.Stop(); sw_elapsed = sw.Elapsed; t2 = System.Convert.ToDouble(sw_elapsed.Ticks); sw.Start(); System.Console.WriteLine(t2 / 10000000);
            // Here is the main subgradient loop
            do
            {
                for (int i = 0; i < N; i++)
                {
                    y_i[i] = 0;
                    for (int r = 0; r < sample_size; r++)
                    {
                        if (1 - Beta_ir[i, r] > 0)      // solve x_ir's with inspection
                        {
                            x_ir[i, r] = 1;
                            temp_UB = temp_UB + 1 - Beta_ir[i, r];
                        }
                        else
                            x_ir[i, r] = 0;
                    }
                }
                //System.Console.WriteLine("from x: " + temp_UB);
                //System.Console.WriteLine("Finished Step-1, calculated x_ir");
                // sw.Stop(); sw_elapsed = sw.Elapsed; t3 = System.Convert.ToDouble(sw_elapsed.Ticks); sw.Start(); System.Console.WriteLine(t3 / 10000000);
                int the_node = 0;
                
                // Solve y's with inspection (Sort according to C_i values and take top-k)
                Dictionary<int, double> c_i_sorted = new Dictionary<int, double>();
                for (int i = 0; i < y_iter.Count; i++)
                {
                    the_node = y_iter[i];
                    c_i_sorted.Add(the_node, c_i[the_node]);
                }

                foreach (KeyValuePair<int, double> item in c_i_sorted.OrderByDescending(x => x.Value).Take(k_init))
                {
                    y_i[item.Key] = 1;
                    y_list.Add(item.Key);
                    temp_UB = temp_UB + c_i[item.Key];
                };

                if (y_list.All(temp_y.Contains))
                {
                    System.Console.Write("DN");  // if we have seen this seed set in the previous iteration then DN -> did not change, no need to calculate the LB again.
                }
                else
                {
                    //find a feasible solution by counting distinct nodes that are accessible (its forward set) with this seed set
                    // Here usually a small number of nodes appear repeatedly, so I keep their Forward set in memory
                    List<HashSet<int>> x_list_int = new List<HashSet<int>>();
                    for (int r = 0; r < sample_size; r++)
                    {
                        x_list_int.Add(new HashSet<int>());
                    }

                    for (int i = 0; i < y_list.Count; i++)
                    {
                        if (y_ever.IndexOf(y_list[i]) > -1)  // if it is in my memory directly add it
                        {
                            int the_index = y_ever.IndexOf(y_list[i]);
                            for (int r = 0; r < sample_size; r++)
                            {
                                foreach (int node in y_influencees[the_index][r])   //my predecessors can activate me
                                {
                                    x_list_int[r].Add(node);
                                }
                            }
                        }
                        else
                        {   // if it is not in my memory do a forward BFS to find all nodes that I can activate
                            y_ever.Add(y_list[i]);
                            int the_index = y_ever.IndexOf(y_list[i]);
                            y_influencees.Add(new List<List<int>>());

                            for (int r = 0; r < sample_size; r++)
                            {
                                y_influencees[the_index].Add(new List<int>());
                                y_influencees[the_index][r].Add(y_list[i]);
                                x_list_int[r].Add(y_list[i]);
                                List<int> single_list = new List<int>();
                                foreach (int node in BFS2(y_list[i], r, 0, SAA_neigh_list).ToList())
                                {
                                    y_influencees[the_index][r].Add(node);
                                    x_list_int[r].Add(node);
                                }
                            }
                        }
                    }
                    temp_LB = temp_LB + x_list_int.Select(x => x.Count).Sum(); // MY lowerbound  or feasible solution value
                   
                }

                if (temp_LB >= z_LB)
                {
                    z_LB = temp_LB;
                    y_best.Clear();
                    for (int jj = 0; jj < y_list.Count; jj++)
                        y_best.Add(y_list[jj]);
                }
               
                if (temp_UB < z_UB)
                {
                    z_UB = temp_UB;
                }

                if (lambda < stopping_threshold || z_LB * (1.0 + stopping_threshold) > z_UB || iter > iterlim)
                {
                    LR_continue = 0;
                }

                // Subgradient Part... First compute all subgradients, then determine step size and update Lagrange multipliers
                double sum_sqr = 0;
           
                for (int i = 0; i < N; i++)
                {
                    for (int r = 0; r < sample_size; r++)
                    {
                        gradient[i, r] = gradient[i, r] + (int)x_ir[i, r];
                        if (y_i[i] == 1)
                        {
                            int the_index = y_ever.IndexOf(i);
                            foreach (int node in y_influencees[the_index][r].ToList())   //my predecessors can activate me
                            {
                                gradient[node, r] = gradient[node, r] - 1;
                            }
                        }
                    }
                }

                for (int i = 0; i < N; i++)
                {
                    for (int r = 0; r < sample_size; r++)
                    {
                        sum_sqr = sum_sqr + gradient[i, r] * gradient[i, r];
                    }
                }
                if (sum_sqr == 0)
                    System.Console.WriteLine("S_SQR=0");
                
                step_size = lambda * (1.05 * z_UB - z_LB) / sum_sqr;
                //System.Console.WriteLine("Finished Step-4, Calculating gradients & stepsize");
                double beta_new = -1;
                int counter = 0;

                List<List<int>> temp_list = new List<List<int>>();

                for (int r = 0; r < sample_size; r++)
                {
                    int[] temp_tree = new int[N];
                    for (int i = 0; i < N; i++)
                    {
                        temp_list.Add(new List<int>());
                        {
                            beta_new = Beta_ir[i, r] + step_size * gradient[i, r];
                            if (beta_new < 0)
                                beta_new = 0;
                            if (beta_new >= 1 && remove_nodes == 0)
                                beta_new = 1;
                            if (beta_new != Beta_ir[i, r])
                            {
                                // If the Lagrange coefficient has changed then update all c_i's that contains this coefficient.
                                // For this we need the RR-list of node i (in all samples)
                                {
                                    c_i[i] = c_i[i] + beta_new - Beta_ir[i, r]; counter++;
                                    List<int> current_list = new List<int>();
                                    if (tree[r, i] == 0)
                                    {
                                        {
                                            temp_list[i].Add(i);
                                            temp_list[i].AddRange(BFS_IM_simple2(i, r, temp_tree, SAA_pred_list, temp_list));
                                            temp_tree[i] = 1;
                                            current_list.AddRange(temp_list[i]);
                                        }
                                    }
                                    else
                                    {
                                        current_list.AddRange(SAA_tree_list[r][i]);
                                    }

                                    foreach (int node in current_list.Skip(1).ToList())
                                    {
                                        c_i[node] = c_i[node] + beta_new - Beta_ir[i, r]; //counter_mem++;
                                    }
                                    counter++; current_list.Clear();
                                }
                                Beta_ir[i, r] = beta_new;
                            }
                        }
                    }
                    temp_list.Clear();
                }

                if (z_UB == z_UB_prev)
                    no_improvement++;

                z_UB_prev = z_UB;
                if (no_improvement >= max_no_improvement_limit)
                {
                    lambda = lambda / 2;
                    no_improvement = 0;
                }

                sw.Stop(); sw_elapsed = sw.Elapsed; t5 = System.Convert.ToDouble(sw_elapsed.Ticks); sw.Start();
                System.Console.WriteLine(iter + "-> " + z_UB + " --- " + z_LB + " --- " + temp_LB + " - " + (z_UB / z_LB - 1) + " - " + t5 / 10000000 + " // y_count:" + y_iter.Count + "( " + counter + ")..." + string.Join(",", y_list.ToArray()));
                if (lambda < stopping_threshold || z_LB * (1.0 + stopping_threshold) > z_UB || iter > iterlim || (t5 / 10000000) > 3600 || (z_UB - z_LB) < 0.00001)
                {
                    LR_continue = 0;
                    //StreamWriter swt = new StreamWriter("betas");
                    //for (int i = 0; i < N; i++)
                    //{
                    //    for (int r = 0; r < sample_size; r++)
                    //    {
                    //        if (SAA_tree_list[r][i].Count > 0)        /// for simplified LR
                    //        {
                    //            if (Beta_ir[i, r] > 0)
                    //                swt.WriteLine("B_" + i + "_" + r + " : " + Beta_ir[i, r]);
                    //        }
                    //    }
                    //}
                    //swt.Close();
                }
                else
                {
                    temp_y = check_list_equality(y_list, temp_y);
                    y_list.Clear();
                    temp_LB = 0;
                    temp_UB = 0;
                    iter++;
                }
                // sw.Stop(); sw_elapsed = sw.Elapsed; t6 = System.Convert.ToDouble(sw_elapsed.Ticks);
                // sw.Start();
                //if (iter < 20)
                //    LR_continue = 1;
                no_of_LR_iter = iter;
            } while (LR_continue > 0);


            solution = z_LB;
            solution_LR = "";
            for (int i = 0; i < y_best.Count; i++)
            {
                solution_LR = solution_LR + ";" + node_set[y_best[i]];
            }
            solution_LR = solution_LR + ";";
            List<int> tempresult = new List<int>();
            for (int j2 = 0; j2 < y_best.Count; j2++)
            {
                tempresult.Add((int)y_best[j2]);
            }
            SAA_LR_list.Add(tempresult);

            z_LRUB = z_UB / sample_size;

            return solution;
        }


        public static double LR_Tiptop(int SAA_m, int remove_nodes, int sample_size)
        {
            double solution = 0.0;
            no_of_LR_iter = 0;
            double[] Beta_r = new double[sample_size];
            double[] x_r = new double[sample_size];
            y_i = new double[N];
            double[] c_i = new double[N];
            double[] p_i = new double[N];
            int j = 0;
            double z_UB = 50.0 * BIGM; //min UB
            double z_LB = -50.0 * BIGM;  // feas sol.
            double lambda = 2;
            int[] gradient = new int[sample_size];
            double temp_UB = 0;
            double temp_LB = 0;
            double step_size = 0;
            double stopping_threshold = 0.00005;
            int max_no_improvement_limit = 30;
            int LR_continue = 1;
            double z_UB_prev = 0;
            int no_improvement = 0;
            int iterlim = 1500;
            int iter = 1;
            int[,] tree = new int[sample_size, N];
            List<int> y_list = new List<int>();
            List<int> y_iter = new List<int>();
            List<int> y_best = new List<int>();
            List<int> y_toremove = new List<int>();
            List<int> temp_y = new List<int>();
            List<int> y_ever = new List<int>();
            List<List<int>> y_no_good = new List<List<int>>();
            List<HashSet<int>> y_influencees = new List<HashSet<int>>();
            double t1, t2, t3, t4, t5, t6;
            sw.Reset();
            if (fixed_probability > 0.1 || N > 10000)
            {
                iterlim = 1000;
                stopping_threshold = 0.001;
                max_no_improvement_limit = 30;
            }
            //sw.Stop(); sw_elapsed = sw.Elapsed; t1 = System.Convert.ToDouble(sw_elapsed.Ticks); sw.Start();
            //initialize beta_ir
            for (int i = 0; i < N; i++)
            {
                y_iter.Add(i);

                c_i[i] = node_score[i] + y_obj[i];
                p_i[i] = N * SAA_N1 * 2;        /// FOR Simplified LR
            }

            for (int r = 0; r < sample_size; r++)
            {
                Beta_r[r] = 1;
            }
            System.Console.WriteLine("Finished initialization of LR! SAA" + SAA_m);

            //sw.Stop(); sw_elapsed = sw.Elapsed; t2 = System.Convert.ToDouble(sw_elapsed.Ticks); sw.Start(); System.Console.WriteLine(t2 / 10000000);

            do
            {
                //int tempctr = 0;
                //stopping_threshold = 0.002;
                for (int i = 0; i < N; i++)
                {
                    //c_i[i] = y_obj[i];        /// FOR Simplified LR
                    y_i[i] = 0;
                }

                for (int r = 0; r < sample_size; r++)
                {
                    //if (x_exist[r, i] == true)    /// FOR simplified LR
                    {
                        if (1 - Beta_r[r] > 0)
                        {
                            x_r[r] = 1;
                            temp_UB = temp_UB + 1 - Beta_r[r];
                            //tempctr++;
                        }
                        else
                            x_r[r] = 0;
                    }
                }
                //System.Console.WriteLine("from x: " + temp_UB);
                //System.Console.WriteLine("Finished Step-1, calculated x_ir");
                // sw.Stop(); sw_elapsed = sw.Elapsed; t3 = System.Convert.ToDouble(sw_elapsed.Ticks); sw.Start(); System.Console.WriteLine(t3 / 10000000);
                int the_node = 0;

                Dictionary<int, double> c_i_sorted = new Dictionary<int, double>();
                for (int i = 0; i < y_iter.Count; i++)
                {
                    the_node = y_iter[i];
                    c_i_sorted.Add(the_node, c_i[the_node]);
                }

                // algorithmic no good cuts
                //if (LR_NG_Cuts == 1)
                //{
                //    bool continueloop = false;
                //    List<int> temp_y_list = new List<int>();
                //    do
                //    {
                //        temp_y_list.Clear();
                //        temp_y_list.AddRange(c_i_sorted.OrderByDescending(x => x.Value).Take(k_init).Select(x => x.Key).ToList());

                //        foreach (List<int> thelist in y_no_good)
                //        {
                //            if (!thelist.Except(temp_y_list).Any())
                //            {
                //                c_i_sorted[temp_y_list[k_init - 1]] = -100000;// * c_i[temp_y_list[k_init - 1]];
                //                continueloop = true;
                //                break;
                //            }
                //            else
                //            {
                //                continueloop = false;
                //            }
                //        }

                //    } while (continueloop);
                //}


                foreach (KeyValuePair<int, double> item in c_i_sorted.OrderByDescending(x => x.Value).Take(k_init))
                {
                    y_i[item.Key] = 1;
                    y_list.Add(item.Key);
                    temp_UB = temp_UB + c_i[item.Key];
                };

                if (y_list.All(temp_y.Contains))
                {
                    System.Console.Write("DN");
                }
                else
                {
                    //find a feasible solution
                    HashSet<int> x_list_int = new HashSet<int>();

                    for (int i = 0; i < y_list.Count; i++)
                    {
                        //temp_LB = temp_LB + y_obj[y_list[i]]; /// for simplified LR
                        if (y_ever.IndexOf(y_list[i]) > -1)
                        {
                            int the_index = y_ever.IndexOf(y_list[i]);

                            foreach (int node in y_influencees[the_index])   //my predecessors can activate me
                            {
                                x_list_int.Add(node);
                            }
                        }
                        else
                        {
                            y_ever.Add(y_list[i]);
                            int the_index = y_ever.IndexOf(y_list[i]);
                            y_influencees.Add(new HashSet<int>());

                            //for (int r = 0; r < sample_size; r++)
                            {
                                //y_influencees[the_index].Add(new List<int>());
                                //y_influencees[the_index].Add(y_list[i]);
                                //x_list_int.Add(y_list[i]);
                                var single_list = hypernodes.Select((x, v) => new { x, v }).Where(a => a.x.IndexOf(y_list[i]) >= 0).Select(a => a.v).ToList();
                                foreach (var node in single_list)
                                {
                                    //if (x_exist[r, node] == true)
                                    {
                                        y_influencees[the_index].Add(Convert.ToUInt16(node));
                                        x_list_int.Add(Convert.ToUInt16(node));
                                    }

                                }
                            }
                        }
                    }
                    temp_LB = temp_LB + x_list_int.Count;
                    // You have the y_obj below for the feasible solution if you need uncomment.


                    for (int i = 0; i < y_list.Count; i++)
                    {
                        int the_index = y_list[i];
                        temp_LB = temp_LB + y_obj[the_index];
                    }
                    //System.Console.WriteLine("Finished Step-3, finding a feasible solution");
                    // sw.Stop(); sw_elapsed = sw.Elapsed; t4 = System.Convert.ToDouble(sw_elapsed.Ticks); sw.Start(); System.Console.WriteLine(t4 / 10000000);
                }

                if (temp_LB >= z_LB)
                {
                    z_LB = temp_LB;
                    y_best.Clear();
                    for (int jj = 0; jj < y_list.Count; jj++)
                        y_best.Add(y_list[jj]);
                }
                else
                {
                    if (LR_NG_Cuts == 1)
                    {
                        y_no_good.Add(new List<int>());
                        System.Console.Write("No good : ");
                        foreach (int node in y_list)
                        {
                            y_no_good[y_no_good.Count - 1].Add(node);
                            System.Console.Write(node + ",");
                        }
                    }
                }

                if (temp_UB < z_UB)
                {
                    z_UB = temp_UB;
                }

                if (lambda < stopping_threshold || z_LB * (1.0 + stopping_threshold) > z_UB || iter > iterlim)
                {
                    LR_continue = 0;
                }
                double sum_sqr = 0;
                // find nodes to delete
                List<int> y_thisround = new List<int>();

                for (int r = 0; r < sample_size; r++)
                {
                    gradient[r] = 0;
                }

                for (int r = 0; r < sample_size; r++)
                {
                    gradient[r] = gradient[r] + (int)x_r[r];
                }

                foreach (int y_elmnt in y_list)
                {
                    int the_index = y_ever.IndexOf(y_elmnt);
                    foreach (int node in y_influencees[the_index].ToList())   //my predecessors can activate me
                    {
                        //if (x_exist[r, node] == true)       //for simplified LR
                        gradient[node] = gradient[node] - 1;
                    }
                }

                for (int r = 0; r < sample_size; r++)
                {
                    //if (x_exist[r, i] == true)
                    sum_sqr = sum_sqr + gradient[r] * gradient[r];
                }
                if (sum_sqr == 0)
                    System.Console.WriteLine("S_SQR=0");
                step_size = lambda * (1.05 * z_UB - z_LB) / sum_sqr;
                //System.Console.WriteLine("Finished Step-4, Calculating gradients & stepsize");
                double beta_new = -1;
                int counter = 0;
                int counter_mem = 0;

                List<List<int>> temp_list = new List<List<int>>();

                for (int r = 0; r < sample_size; r++)
                {
                    int[] temp_tree = new int[N];
                    counter_mem = 0;//System.Console.Write(r + "..");

                    temp_list.Add(new List<int>());
                    //if (x_exist[r, i] == true)        /// for simplified LR
                    {
                        beta_new = Beta_r[r] + step_size * gradient[r];
                        if (beta_new < 0)
                            beta_new = 0;
                        if (beta_new >= 1 && remove_nodes == 0)
                            beta_new = 1;
                        if (beta_new != Beta_r[r])
                        {
                            foreach (int node in hypernodes[r].ToList())
                            {
                                c_i[node] = c_i[node] + beta_new - Beta_r[r]; //counter_mem++;
                            }
                            counter++;

                            Beta_r[r] = beta_new;
                        }
                    }

                    temp_list.Clear();
                }

                if (z_UB == z_UB_prev)
                    no_improvement++;

                z_UB_prev = z_UB;
                if (no_improvement >= max_no_improvement_limit)
                {
                    lambda = lambda / 2;
                    no_improvement = 0;
                }

                sw.Stop(); sw_elapsed = sw.Elapsed; t5 = System.Convert.ToDouble(sw_elapsed.Ticks); sw.Start();
                System.Console.WriteLine(iter + "-> " + z_UB + " --- " + z_LB + " --- " + temp_LB + " - " + (z_UB / z_LB - 1) + " - " + t5 / 10000000 + " // y_count:" + y_iter.Count + "( " + counter + ")..." + string.Join(",", y_list.ToArray()));
                if (lambda < stopping_threshold || z_LB * (1.0 + stopping_threshold) > z_UB || iter > iterlim || (t5 / 10000000) > 3600 || (z_UB - z_LB) < 0.00001)
                {
                    LR_continue = 0;
                    //StreamWriter swt = new StreamWriter("betas");
                    //for (int i = 0; i < N; i++)
                    //{
                    //    for (int r = 0; r < sample_size; r++)
                    //    {
                    //        if (SAA_tree_list[r][i].Count > 0)        /// for simplified LR
                    //        {
                    //            if (Beta_ir[i, r] > 0)
                    //                swt.WriteLine("B_" + i + "_" + r + " : " + Beta_ir[i, r]);
                    //        }
                    //    }
                    //}
                    //swt.Close();
                }
                else
                {
                    temp_y = check_list_equality(y_list, temp_y);
                    y_list.Clear();
                    temp_LB = 0;
                    temp_UB = 0;
                    iter++;
                }
                // sw.Stop(); sw_elapsed = sw.Elapsed; t6 = System.Convert.ToDouble(sw_elapsed.Ticks);
                // sw.Start();
                //if (iter < 20)
                //    LR_continue = 1;
                no_of_LR_iter = iter;
            } while (LR_continue > 0);


            solution = z_LB;
            solution_LR = "";
            for (int i = 0; i < y_best.Count; i++)
            {
                solution_LR = solution_LR + ";" + node_set[y_best[i]];
            }
            solution_LR = solution_LR + ";";
            List<int> tempresult = new List<int>();
            for (int j2 = 0; j2 < y_best.Count; j2++)
            {
                tempresult.Add((int)y_best[j2]);
            }
            SAA_LR_list.Add(tempresult);

            z_LRUB = z_UB;

            return solution;
        }
  

        public static List<int> check_list_equality(List<int> y_list, List<int> temp_y)
        {
            var a = y_list.All(temp_y.Contains) && y_list.Count == temp_y.Count;
            if (a != true)
            {
                temp_y.Clear();
                temp_y.AddRange(y_list);
            }
            return temp_y;
        }


        public static List<int> return_UInt16_list(string str)
        {
            List<int> mylist = new List<int>();
            string[] str_arr = str.Split(new char[] { ',' });
            int item;
            for (int i = 0; i < str_arr.Length; i++)
            {
                if (str_arr[i] != "")
                {
                    item = Convert.ToUInt16(str_arr[i]);
                    mylist.Add(item);
                }
            }
            return mylist;
        }

        public static double evaluate_influence_fast(List<int> initial_set, int SAA_m, int SAA_n2, int type) //evaluates from nodeIDs not actual nodes!!!
        {
            double influence = -1;

            return influence;
        }

        public static double evaluate_influence(List<int> initial_set, int SAA_m, int SAA_n2, int type) //evaluates from nodeIDs not actual nodes!!!
        {
            int temp = 0;
            double influence = 0;
            counter_prob = 0;
            counter_prob2 = 0;
            //only_arcs = unique_arcs.Select(b =>
            //               new ArcsUInt16Weighted { Head = b.Key.Head, Tail = b.Key.Tail, Weight = b.W }).AsParallel().ToList();
            //var left_nodes_neighbours = unique_arcs.GroupBy(a => a.Key.Head).Select(a => new { Head = node_index[(int)a.Key], RefList = a.Select(b => node_index[(int)b.Key.Tail]).ToList() }).AsParallel().ToList();
            var g1 = unique_arcs.GroupBy(a => a.Key.Head).Select(a => new { TailList = a.Select(b => node_index[(int)b.Key.Tail]).ToList(), WList = a.Select(b => b.W).ToList() }).AsParallel().ToList();
            List<List<int>> arclist = new List<List<int>>();
            List<double> wlist = new List<double>();
            List<int> tails = new List<int>();
            int list_size = unique_arcs.Count;
            int head = 0;
            int tail = 0;
            Random myrand = new Random();


            for (int i = 0; i < N; i++)
                arclist.Add(new List<int>());

            for (int i = 0; i < list_size; i++)
            {
                head = node_index[(int)unique_arcs[i].Key.Head];
                tail = node_index[(int)unique_arcs[i].Key.Tail];
                arclist[head].Add(i);
                tails.Add(tail);
                wlist.Add(unique_arcs[i].W);
            }

            int temp_count = 0;
            //swtemp2 = new StreamWriter("prob2.txt");
            for (int t = 0; t < SAA_n2; t++)
            {
                if (t % 1000 == 0)
                {
                    System.Console.Write("-" + t);
                }
                if (type == 3)
                {
                    independent_cascade_sub_for_greedy(initial_set, SAA_m, t, type);
                    temp = temp + active_set.Count;
                }

                else
                //independent_cascade_sub(initial_set, SAA_m, t, type);
                {
                    temp = temp + independent_cascade_sub_fast2(initial_set, SAA_m, t, type, arclist, wlist, tails, myrand);
                    //int temp2= independent_cascade_sub_fast2(initial_set, SAA_m, t, type, arclist, wlist, tails, myrand);
                }

            }
            influence = (double)temp / SAA_n2;
            System.Console.WriteLine();
            System.Console.WriteLine("cp:" + counter_prob + "/" + counter_prob2);
            //swtemp2.Close();
            return influence;
        }


        public static int independent_cascade_sub_fast(List<int> IS, int SAA_m, int t, int type, List<List<int>> arclist, List<double> wlist, List<int> tails, Random myrand)
        {
            var vertices = new Queue<int>();
            int tail = 0;
            double prop;
            bool[] marked = new bool[N];

            foreach (int element in IS.ToList())
                vertices.Enqueue(element);

            HashSet<int> visited_nodes = new HashSet<int>();
            //            bool[] marked = new bool[N];
            while (vertices.Any())
            {
                int curVertex = vertices.Dequeue();
                visited_nodes.Add(curVertex);
                marked[curVertex] = true;
                foreach (int arc in arclist[curVertex].ToList())
                {
                    tail = tails[arc];
                    if (marked[tail] == false)
                    {
                        prop = myrand.NextDouble();
                        if (prop <= wlist[arc])
                        {
                            marked[tail] = true;
                            visited_nodes.Add(tail);
                            vertices.Enqueue(tail);
                        }
                    }
                }
            }
            return visited_nodes.Count;
        }

        public static int independent_cascade_sub_fast2(List<int> IS, int SAA_m, int t, int type, List<List<int>> arclist, List<double> wlist, List<int> tails, Random myrand)
        {
            active_set = new List<int>(IS);
            CA_set = new List<int>(IS);
            NA_set = new List<int>();
            bool[] marked = new bool[N];

            foreach (int element in IS.ToList())
                marked[element] = true;
            int curr_neigh;
            int temp_neigh_index = -1;
            double temp_arc_weight = 0;
            do
            {
                temp_neigh_index = CA_set[0];
                //swvar.WriteLine(CA_set[i]+" ("+t+")");
                foreach (int edge in arclist[temp_neigh_index].ToList())   //tüm komşuları aktifleştirmeye çalış
                {
                    curr_neigh = tails[edge];
                    if (!marked[curr_neigh]) //tabi neighbour active set içinde değilse aktifleştirmeye çalış
                    {
                        temp_arc_weight = myrand.NextDouble();
                        //swtemp2.WriteLine(temp_arc_weight);
                        if (temp_arc_weight <= wlist[edge])
                        {
                            NA_set.Add(curr_neigh);
                            marked[curr_neigh] = true;
                            counter_prob++;
                            //swvar3.WriteLine(curr_neigh + " <- act by " + temp_neigh_index + " (" + t + ") prb:" + temp_arc_weight);
                        }
                    }
                }
                CA_set.RemoveAt(0);
                active_set.AddRange(NA_set);
                CA_set.AddRange(NA_set);
                NA_set.Clear();
            } while (CA_set.Count > 0);
            return active_set.Count;
        }

        public static void independent_cascade_sub(List<int> IS, int SAA_m, int t, int type)
        {
            string[] neigh_arr;
            active_set = new List<int>(IS);
            CA_set = new List<int>(IS);
            NA_set = new List<int>();
            int curr_neigh;
            int curr_neigh_index;
            string temp_arc;

            int temp_arc_index;
            int temp_neigh_index = -1;
            int arcID = -1;
            double temp_arc_weight = 0;
            do
            {
                //foreach (int i in CA_set.ToList())
                int i = 0;
                {
                    temp_neigh_index = System.Convert.ToUInt16(CA_set[i]);
                    //swvar.WriteLine(CA_set[i]+" ("+t+")");
                    if (neigh_edge_list[temp_neigh_index] != null)
                    {
                        foreach (int edge in neigh_edge_list[temp_neigh_index])   //tüm komşuları aktifleştirmeye çalış
                        {
                            //arcID = (int)edge;
                            curr_neigh = (int)arcs_intID[edge, 1];
                            if (!active_set.Contains(curr_neigh)) //tabi neighbour active set içinde değilse aktifleştirmeye çalış
                            {
                                //temp_arc = CA_set[i].ToString() + "," + curr_neigh.ToString();
                                //if (influence_prob[CA_set[i], curr_neigh] >= prop_prob)
                                //temp_arc_index = arcID_set.IndexOf(temp_arc);
                                if (type == 1)
                                    temp_arc_weight = GetUniform();// SAA_2_prb[arcID, t];
                                //if (type == 0)
                                //    temp_arc_weight = SAA_prob[arcID, t, SAA_m];
                                if (type == 2)
                                    temp_arc_weight = GetUniform(); //SAA_2_prob[arcID, t];

                                if (type == 3)
                                    temp_arc_weight = SAA_3_prob[edge, t];

                                if (type == 4)  // for only evaluate
                                    temp_arc_weight = SAA_3_prob[edge, t];

                                if (type == 5)  // for only evaluate
                                    temp_arc_weight = SAA_3_prob_dbl[edge][t];

                                if (type == 6)  // for only evaluate
                                    temp_arc_weight = SAA_float[edge, t];

                                //swtemp2.WriteLine(temp_arc_weight);
                                counter_prob2++;
                                prop_prob = 1 - weight_set[edge];
                                if (temp_arc_weight >= prop_prob)
                                {
                                    NA_set.Add(curr_neigh);
                                    counter_prob++;
                                    //swvar3.WriteLine(curr_neigh + " <- act by " + temp_neigh_index + " (" + t + ") prb:" + temp_arc_weight);
                                }
                            }
                        }
                    }
                    CA_set.Remove(CA_set[i]);
                    foreach (int j in NA_set)
                    {
                        active_set.Add(j);
                        CA_set.Add(j);
                    }
                    NA_set.Clear();
                    i++;
                }
            } while (CA_set.Count > 0);
        }

        public static void independent_cascade_sub_for_greedy(List<int> IS, int SAA_m, int t, int type)
        {
            string[] neigh_arr;
            active_set = new List<int>(IS);
            CA_set = new List<int>(IS);
            NA_set = new List<int>();
            int curr_neigh;

            //SAA_neigh = new string[N, sample_size, SAA_M];
            //SAA_neigh[temp_neigh_index,t,SAA_m]

            int temp_neigh_index = -1;
            int arcID = -1;
            double temp_arc_weight = 0;
            do
            {
                //foreach (int i in CA_set.ToList())
                int i = 0;
                {
                    temp_neigh_index = System.Convert.ToUInt16(CA_set[i]);
                    //swvar.WriteLine(CA_set[i]+" ("+t+")");

                    if (SAA_neigh[temp_neigh_index, t, SAA_m] != null)
                    {
                        neigh_arr = SAA_neigh[temp_neigh_index, t, SAA_m].Split(new char[] { ',' });     //neighbour stringini split edelim
                        foreach (string str in neigh_arr)   //tüm komşuları aktifleştirmeye çalış
                        {
                            curr_neigh = Convert.ToUInt16(str);
                            if (!active_set.Contains(curr_neigh)) //tabi neighbour active set içinde değilse aktifleştirmeye çalış
                            {
                                NA_set.Add(curr_neigh);
                            }
                        }
                    }
                    CA_set.Remove(CA_set[i]);
                    foreach (int j in NA_set)
                    {
                        active_set.Add(j);
                        CA_set.Add(j);
                    }
                    NA_set.Clear();
                    i++;
                }
            } while (CA_set.Count > 0);
        }



        #region statistics

        // The random generator seed can be set three ways:
        // 1) specifying two non-zero unsigned integers
        // 2) specifying one non-zero unsigned integer and taking a default value for the second
        // 3) setting the seed from the system time

        public static void SetSeed(uint u, uint v)
        {

            if (u != 0) m_w = u;
            if (v != 0) m_z = v;
        }

        public static void SetSeed(uint u)
        {
            m_w = u;
        }

        public static void SetSeedFromSystemTime()
        {
            System.DateTime dt = System.DateTime.Now;
            long x = dt.ToFileTime();
            SetSeed((uint)(x >> 16), (uint)(x % 4294967296));
        }

        // Produce a uniform random sample from the open interval (0, 1).
        // The method will not return either end point.
        public static double GetUniform()
        {
            // 0 <= u < 2^32
            uint u = GetUint();
            // The magic number below is 1/(2^32 + 2).
            // The result is strictly between 0 and 1.
            return (u + 1.0) * 2.328306435454494e-10;
        }

        // This is the heart of the generator.
        // It uses George Marsaglia's MWC algorithm to produce an unsigned integer.
        // See http://www.bobwheeler.com/statistics/Password/MarsagliaPost.txt
        private static uint GetUint()
        {
            m_z = 36969 * (m_z & 65535) + (m_z >> 16);
            m_w = 18000 * (m_w & 65535) + (m_w >> 16);
            return (m_z << 16) + m_w;
        }

        // Get normal (Gaussian) random sample with mean 0 and standard deviation 1
        public static double GetNormal()
        {
            // Use Box-Muller algorithm
            double u1 = GetUniform();
            double u2 = GetUniform();
            double r = Math.Sqrt(-2.0 * Math.Log(u1));
            double theta = 2.0 * Math.PI * u2;
            return r * Math.Sin(theta);
        }

        // Get normal (Gaussian) random sample with specified mean and standard deviation
        public static double GetNormal(double mean, double standardDeviation)
        {
            if (standardDeviation <= 0.0)
            {
                string msg = string.Format("Shape must be positive. Received {0}.", standardDeviation);
                throw new ArgumentOutOfRangeException(msg);
            }
            return mean + standardDeviation * GetNormal();
        }

        // Get exponential random sample with mean 1
        public static double GetExponential()
        {
            return -Math.Log(GetUniform());
        }

        // Get exponential random sample with specified mean
        public static double GetExponential(double mean)
        {
            if (mean <= 0.0)
            {
                string msg = string.Format("Mean must be positive. Received {0}.", mean);
                throw new ArgumentOutOfRangeException(msg);
            }
            return mean * GetExponential();
        }

        public static double GetGamma(double shape, double scale)
        {
            // Implementation based on "A Simple Method for Generating Gamma Variables"
            // by George Marsaglia and Wai Wan Tsang.  ACM Transactions on Mathematical Software
            // Vol 26, No 3, September 2000, pages 363-372.

            double d, c, x, xsquared, v, u;

            if (shape >= 1.0)
            {
                d = shape - 1.0 / 3.0;
                c = 1.0 / Math.Sqrt(9.0 * d);
                for (; ; )
                {
                    do
                    {
                        x = GetNormal();
                        v = 1.0 + c * x;
                    }
                    while (v <= 0.0);
                    v = v * v * v;
                    u = GetUniform();
                    xsquared = x * x;
                    if (u < 1.0 - .0331 * xsquared * xsquared || Math.Log(u) < 0.5 * xsquared + d * (1.0 - v + Math.Log(v)))
                        return scale * d * v;
                }
            }
            else if (shape <= 0.0)
            {
                string msg = string.Format("Shape must be positive. Received {0}.", shape);
                throw new ArgumentOutOfRangeException(msg);
            }
            else
            {
                double g = GetGamma(shape + 1.0, 1.0);
                double w = GetUniform();
                return scale * g * Math.Pow(w, 1.0 / shape);
            }
        }

        public static double GetChiSquare(double degreesOfFreedom)
        {
            // A chi squared distribution with n degrees of freedom
            // is a gamma distribution with shape n/2 and scale 2.
            return GetGamma(0.5 * degreesOfFreedom, 2.0);
        }

        public static double GetInverseGamma(double shape, double scale)
        {
            // If X is gamma(shape, scale) then
            // 1/Y is inverse gamma(shape, 1/scale)
            return 1.0 / GetGamma(shape, 1.0 / scale);
        }

        public static double GetWeibull(double shape, double scale)
        {
            if (shape <= 0.0 || scale <= 0.0)
            {
                string msg = string.Format("Shape and scale parameters must be positive. Recieved shape {0} and scale{1}.", shape, scale);
                throw new ArgumentOutOfRangeException(msg);
            }
            return scale * Math.Pow(-Math.Log(GetUniform()), 1.0 / shape);
        }

        public static double GetCauchy(double median, double scale)
        {
            if (scale <= 0)
            {
                string msg = string.Format("Scale must be positive. Received {0}.", scale);
                throw new ArgumentException(msg);
            }

            double p = GetUniform();

            // Apply inverse of the Cauchy distribution function to a uniform
            return median + scale * Math.Tan(Math.PI * (p - 0.5));
        }

        public static double GetStudentT(double degreesOfFreedom)
        {
            if (degreesOfFreedom <= 0)
            {
                string msg = string.Format("Degrees of freedom must be positive. Received {0}.", degreesOfFreedom);
                throw new ArgumentException(msg);
            }

            // See Seminumerical Algorithms by Knuth
            double y1 = GetNormal();
            double y2 = GetChiSquare(degreesOfFreedom);
            return y1 / Math.Sqrt(y2 / degreesOfFreedom);
        }

        // The Laplace distribution is also known as the double exponential distribution.
        public static double GetLaplace(double mean, double scale)
        {
            double u = GetUniform();
            return (u < 0.5) ?
                mean + scale * Math.Log(2.0 * u) :
                mean - scale * Math.Log(2 * (1 - u));
        }

        public static double GetLogNormal(double mu, double sigma)
        {
            return Math.Exp(GetNormal(mu, sigma));
        }

        public static double GetBeta(double a, double b)
        {
            if (a <= 0.0 || b <= 0.0)
            {
                string msg = string.Format("Beta parameters must be positive. Received {0} and {1}.", a, b);
                throw new ArgumentOutOfRangeException(msg);
            }

            // There are more efficient methods for generating beta samples.
            // However such methods are a little more efficient and much more complicated.
            // For an explanation of why the following method works, see
            // http://www.johndcook.com/distribution_chart.html#gamma_beta

            double u = GetGamma(a, 1.0);
            double v = GetGamma(b, 1.0);
            return u / (u + v);
        }
        #endregion

    }
}
