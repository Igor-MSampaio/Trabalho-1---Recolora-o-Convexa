using Gurobi;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;

namespace RecoloracaoConvexa
{
    class Program
    {
        public class Input
        {
            public int nVertices;
            public int nCores;
            public List<int> caminhoColorido;
        }

        public class BB_Node
        {
            public Tuple<int, int> varToFix;
            public bool fix1;
            public BB_Node parent;
            public BB_Node left_child;
            public BB_Node right_child;
            public double UB;

            public BB_Node(Tuple<int, int> varToFix, bool fix1, BB_Node parent, double UB)
            {
                this.varToFix = varToFix;
                this.fix1 = fix1;
                this.parent = parent;
                this.left_child = null;
                this.right_child = null;
                this.UB = UB;
            }
        }

        public class Config
        {
            public bool turn_heuristic_on;
            public bool turn_cuts_on;
            public int num_max_cuts_per_iter;
            public double violation_threshold;
            public double time_limit;

            public Config()
            {
                this.turn_heuristic_on = true;
                this.turn_cuts_on = true;
                this.num_max_cuts_per_iter = 10;
                this.violation_threshold = 0.001;
                this.time_limit = 1800;
            }
        }

        public class Statistics
        {
            public int n_BB_nodes;
            public double heuristic_sol_val;
            public double objective_value;
            public double first_linear_relax_val;
            public double spent_time_in_seconds;
            public int[,] solucao_heuristica;
 
            public Statistics()
            {
                this.n_BB_nodes = 0;
                this.heuristic_sol_val = 0;
                this.objective_value = 0;
                this.first_linear_relax_val = 0;
                this.spent_time_in_seconds = 0;
            }
        }

        public static double e = 0.000001;
        public static  Input input = new Input();
        public static GRBEnv env;
        public static  GRBModel model;
        public static  GRBVar[,] X;
        public static Config config = new Config();
        public static Statistics statistics = new Statistics();

        public static void read_input(string instancia)
        {
            StreamReader file = File.OpenText(@"Instancias\" + instancia + ".txt");
            string[] words = file.ReadLine().Split(' ');
            input.nVertices = int.Parse(words[0]);
            input.nCores = int.Parse(words[1]);

            input.caminhoColorido = new List<int>();

            string line;

            while ((line = file.ReadLine()) != null)
            {
                input.caminhoColorido.Add(int.Parse(line));
            }
        }

        public static void build_initial_model(bool lp_relax)
        {
            env = new GRBEnv(true);
            env.Set("LogFile", "Recoloracao.log");
            env.Set(GRB.IntParam.LogToConsole, 0);
            env.Start();

            model = new GRBModel(env);
            model.ModelName = "Recoloracao";
            model.Set(GRB.IntParam.Threads, 1);
            model.Set(GRB.IntParam.Method, 1);
            model.Set(GRB.DoubleParam.TimeLimit, 1800);

            X = new GRBVar[input.nVertices, input.nCores];
            var obj = new GRBLinExpr();

            for (int i = 0; i < input.nVertices; i++)
            {
                for (int j = 0; j < input.nCores; j++)
                {
                    string nome = "X_ " + i.ToString() + "_c " + j.ToString();

                    if (lp_relax)
                        X[i, j] = model.AddVar(0.0, 1.0, 1.0, GRB.CONTINUOUS, nome);
                    else
                        X[i, j] = model.AddVar(0.0, 1.0, 1.0, GRB.BINARY, nome);

                    if (input.caminhoColorido[i] != j + 1)
                        obj.AddTerm(1, X[i, j]);
                }
            }

            model.SetObjective(obj, GRB.MINIMIZE);

            // Restrição 1
            for (int i = 0; i < input.nVertices; i++)
            {
                GRBLinExpr soma = new GRBLinExpr();

                for (int c = 0; c < input.nCores; c++)
                    soma.AddTerm(1.0, X[i, c]);

                model.AddConstr(soma, GRB.EQUAL, 1, "Cor vértice" + i);
            }

            // Restrição 2
            for (int p = 0; p < input.nVertices - 2; p++)
            {
                for (int r = p + 2; r < input.nVertices; r++)
                {
                    for (int q = p + 1; q < r; q++)
                    {
                        for (int k = 0; k < input.nCores; k++)
                        {
                            model.AddConstr(X[p, k] - X[q, k] + X[r, k], GRB.LESS_EQUAL, 1, "p" + p + "-q" + q + "-r" + r + "-k" + k);
                        }
                    }
                }
            }
        }

        public static void Heuristica()
        {
            List<int> contCores = new List<int>();

            for (int i = 0; i < input.nCores; i++)
            {
                contCores.Add(0);
            }

            foreach (var item in input.caminhoColorido)
            {
                contCores[item - 1]++;
            }

            int corMaior = 0;
            int contMaior = 0;

            for (int i = 0; i < contCores.Count; i++)
            {
                if (contMaior < contCores[i])
                {
                    contMaior = contCores[i];
                    corMaior = i;
                }
                else if (contMaior == contCores[i])
                {
                    if (input.caminhoColorido[0] == corMaior + 1 && input.caminhoColorido[input.caminhoColorido.Count - 1] == corMaior + 1)
                    {
                        contMaior = contCores[i];
                        corMaior = i;
                    }
                    else if (input.caminhoColorido[0] == corMaior + 1 || input.caminhoColorido[input.caminhoColorido.Count - 1] == corMaior + 1)
                    {
                        if (input.caminhoColorido[0] == i + 1 || input.caminhoColorido[input.caminhoColorido.Count - 1] == i + 1)
                        {
                            continue;
                        }
                        else
                        {
                            contMaior = contCores[i];
                            corMaior = i;
                        }
                    }
                }
            }

            int naoTrocaIni = 0;
            int naoTrocaFim = 0;
            int corIni = input.caminhoColorido[0];
            int corFim = input.caminhoColorido[input.caminhoColorido.Count - 1];

            foreach (var item in input.caminhoColorido)
            {
                if ((item - 1) == corMaior || item != corIni)
                    break;
                else
                    naoTrocaIni++;
            }

            for (int i = input.caminhoColorido.Count - 1; i >= 0; i--)
            {
                if ((input.caminhoColorido[i] - 1) == corMaior || input.caminhoColorido[i] != corFim || corIni == corFim)
                    break;
                else
                    naoTrocaFim++;
            }

            List<int> novoCaminho = new List<int>();

            for (int i = 0; i < input.caminhoColorido.Count; i++)
            {
                if (i < naoTrocaIni || input.caminhoColorido.Count - i <= naoTrocaFim)
                    novoCaminho.Add(input.caminhoColorido[i]);
                else
                {
                    novoCaminho.Add(corMaior + 1);
                }
            }

            statistics.solucao_heuristica = new int[input.nVertices, input.nCores];

            for (int i = 0; i < input.nVertices; i++)
            {
                for (int j = 0; j < input.nCores; j++)
                {
                    if (j + 1 == novoCaminho[i])
                        statistics.solucao_heuristica[i, j] = 1;
                    else
                        statistics.solucao_heuristica[i, j] = 0;
                }
            }

            statistics.heuristic_sol_val = (input.caminhoColorido.Count - contMaior - (naoTrocaIni + naoTrocaFim));
        }

        public static int solve_LP()
        {
            model.Optimize();
            return model.Status;
        }

        public static void solve_IP()
        {
            build_initial_model(false);
            model.Optimize();

            Console.WriteLine("--------------------------------");
            Console.WriteLine($"Tempo de execução para a otimização mais recente: {model.Runtime}");
            Console.WriteLine($"Número de nós de ramificação e corte explorados na otimização mais recente: {model.NodeCount}");
            Console.WriteLine($"Valor objetivo para a solução atual: {model.ObjVal}");
            Console.WriteLine($"Número de variáveis = {model.NumVars}");
            Console.WriteLine($"Número de restrições lineares = {model.NumConstrs}");
            Console.WriteLine($"Número de coeficientes diferentes de zero na matriz de restrição = {model.NumNZs}");
            Console.WriteLine("--------------------------------");
        }

        public static void add_BB_node(List<BB_Node> L, BB_Node no)
        {
            L.Add(no);
        }

        public static BB_Node select_BB_node(List<BB_Node> L)
        {
            BB_Node selected_node = null;
            double max_UB = 0.0;

            foreach (BB_Node node in L)
            {
                if (node.UB > max_UB)
                {
                    max_UB = node.UB;
                    selected_node = node;
                }
            }

            L.Remove(selected_node);

            return selected_node;
        }

        public static BB_Node create_BB_node(Tuple<int, int> varToFix, bool fix1, BB_Node parent, double UB)
        {
            statistics.n_BB_nodes += 1;
            BB_Node no = new BB_Node(varToFix, fix1, parent, UB);
            return no;
        }

        public static void clear_branch_contraints()
        {
            for (int i = 0; i < input.nVertices; i++)
            {
                for (int j = 0; j < input.nCores; j++)
                {
                    X[i, j].UB = 1;
                    X[i, j].LB = 0;
                }
            }
        }

        public static void add_branch_contraints(BB_Node n_i)
        {
            if (n_i.parent != null)
            {
                if (n_i.fix1)
                    X[n_i.varToFix.Item1, n_i.varToFix.Item2].LB = 1;
                else
                    X[n_i.varToFix.Item1, n_i.varToFix.Item2].UB = 0;

                add_branch_contraints(n_i.parent);
            }
        }

        public static Tuple<int, int> select_frac_var()
        {
            Tuple<int, int> mais_frac_i_j = new Tuple<int, int>(0, 0);
            double mais_frac = X[0, 0].Get(GRB.DoubleAttr.X);

            for (int i = 0; i < input.nVertices; i++)
            {
                for (int j = 0; j < input.nCores; j++)
                {
                    if (Math.Abs(X[i, j].Get(GRB.DoubleAttr.X) - 0.5) < Math.Abs(mais_frac - 0.5))
                    {
                        mais_frac_i_j = new Tuple<int, int>(i, j);
                        mais_frac = X[i, j].Get(GRB.DoubleAttr.X);
                    }
                }
            }

            if (mais_frac < e || mais_frac > 1 - e)
                return null;
            else
                return mais_frac_i_j;
        }

        public static bool is_integer()
        {
            return select_frac_var() == null;
        }

        public static void copy_solution(int[,] x_int)
        {
            for (int i = 0; i < input.nVertices; i++)
            {
                for (int j = 0; j < input.nCores; j++)
                {
                    if (X[i, j].Get(GRB.DoubleAttr.X) > e)
                        x_int[i, j] = 1;
                    else
                        x_int[i, j] = 0;
                }
            }
        }

        public static double update_UB(BB_Node node)
        {
            if (node.left_child == null || node.right_child == null)
                return node.UB;

            return node.UB = Math.Min(update_UB(node.left_child), update_UB(node.right_child));
        }

        public static void BranchAndBound()
        {
            var stopwatch = new Stopwatch();
            stopwatch.Start();

            build_initial_model(true);
            int[,] best_sol = statistics.solucao_heuristica;
            double z_ = statistics.heuristic_sol_val;

            BB_Node n_0 = create_BB_node(new Tuple<int, int>(-1, -1), false, null, int.MaxValue);
            List<BB_Node> L = new List<BB_Node>();
            add_BB_node(L, n_0);
            var is_root_node = true;
            BB_Node root_node = null;

            while (L.Count != 0)
            {
                
                if (is_root_node == false)
                {
                    update_UB(root_node);

                    if (root_node.UB >= z_ + e)
                    {
                        break;
                    }
                }


                stopwatch.Stop();
                statistics.spent_time_in_seconds = stopwatch.ElapsedMilliseconds;

                if (statistics.spent_time_in_seconds / 1000 > config.time_limit)
                {
                    time_limit_exceeded(best_sol, z_, root_node.UB);
                    return;
                }

                stopwatch.Start();

                BB_Node n_i = select_BB_node(L);
                clear_branch_contraints();
                add_branch_contraints(n_i);

                var status = solve_LP();
                double z_i;

                if (status == 3)
                    continue;
                else
                {
                    z_i = model.ObjVal;
                    n_i.UB = z_i;

                    if (is_root_node)
                    {
                        statistics.first_linear_relax_val = z_i;
                        root_node = n_i;
                        is_root_node = false;
                    }
                }

                if (z_i >= z_ + e)
                    continue;

                if (is_integer())
                {
                    z_ = z_i;
                    copy_solution(best_sol);
                    continue;
                }

                Tuple<int, int> i_var = select_frac_var();
                n_i.left_child = create_BB_node(i_var, true, n_i, z_i);
                n_i.right_child = create_BB_node(i_var, false, n_i, z_i);
                add_BB_node(L, n_i.left_child);
                add_BB_node(L, n_i.right_child);
            }

            stopwatch.Stop();
            statistics.spent_time_in_seconds = stopwatch.Elapsed.TotalSeconds;
            statistics.objective_value = z_;
            print_solution(best_sol);
        }

        public static void time_limit_exceeded(int[,] sol, double z_, double UB)
        {

            Console.WriteLine("--------------------------------");
            Console.WriteLine($"TimeLimit de 30 minutos excedido");
            Console.WriteLine($"Número de nós explorados até o TimeLimit: {statistics.n_BB_nodes}");
            Console.WriteLine($"Valor da relaxação inicial = {statistics.first_linear_relax_val}");
            Console.WriteLine($"Melhor Valor objetivo até o TimeLimit: {statistics.objective_value}");
            Console.WriteLine($"Valor encontrado pela heuristica: {statistics.heuristic_sol_val}");
            Console.WriteLine($"Número de variáveis = {model.NumVars}");
            Console.WriteLine($"Número de restrições lineares = {model.NumConstrs}");
            Console.WriteLine($"Número de coeficientes diferentes de zero na matriz de restrição = {model.NumNZs}");
            Console.WriteLine("--------------------------------");

        }

        public static void print_solution(int[,] sol)
        {
            Console.WriteLine("--------------------------------");
            Console.WriteLine($"Tempo total: {statistics.spent_time_in_seconds}");
            Console.WriteLine($"Número de nós explorados: {statistics.n_BB_nodes}");
            Console.WriteLine($"Valor da relaxação inicial = {statistics.first_linear_relax_val}");
            Console.WriteLine($"Valor objetivo: {statistics.objective_value}");
            Console.WriteLine($"Valor encontrado pela heuristica: {statistics.heuristic_sol_val}");
            Console.WriteLine($"Número de variáveis = {model.NumVars}");
            Console.WriteLine($"Número de restrições lineares = {model.NumConstrs}");
            Console.WriteLine($"Número de coeficientes diferentes de zero na matriz de restrição = {model.NumNZs}");
            Console.WriteLine("--------------------------------");
        }

        public static void Main(string[] args)
        {
            for (int i = 10; i <= 50; i += 10)
            {
                for (int j = 2; j <= 10; j++)
                {
                    var instancia = $"rand_{i}_{j}";
                    read_input(instancia);
                    Console.WriteLine("-------------------------------------------------------------------");
                    Console.WriteLine(instancia);
                    Console.WriteLine("--------------------------------");
                    Console.WriteLine("Resultado Branch And Bound: ");
                    Heuristica();
                    BranchAndBound();
                    Console.WriteLine("Resultado Modelo Inteiro: ");
                    solve_IP();
                }
            }
        }
    }
}
