// g++-9 omp_static_sched_rng_per_thread.cpp -fopenmp -o gcc_omp_static_sched_rng_per_thread
// clang++ -fomit-frame-pointer -Xpreprocessor -fopenmp -m64 -std=c++11 -L/usr/local/opt/libomp/lib -lomp omp_static_sched_rng_per_thread.cpp -o omp_static_sched_rng_per_thread

#include <vector>
#include <random>
#include <iostream>
#include <cstdlib>

#include <omp.h>

std::random_device rd;
std::vector<std::mt19937> mt19937_vec;

// std::vector<double> rval;
double rval[20000];

double UniformRandom(std::mt19937 &mt_gen)
{
	return std::generate_canonical<double, 10>(mt_gen);
}

int main(int argc, char **argv)
{
  int num_threads = 0;
  int tid;
  int num_cells = 20000;
  double rnum;

  std::cout << "argc = " << argc << std::endl;
  // for (int i = 0; i < argc; ++i) 
        // std::cout << argv[i] << std::endl;
  if (argc < 2)
  {
    std::cout << "Provide # threads as an argument." << std::endl;
    exit(-1);
  }

  // std::cout << "omp_get_num_threads() (serial section) = " << omp_get_num_threads() << std::endl;  // = 1 !!
  std::cout << "omp_get_max_threads() (serial section) = " << omp_get_max_threads() << std::endl;
  if (argc > 1)
  {
    num_threads = atoi(argv[1]);
    omp_set_num_threads(num_threads);
  }
  std::cout << "num_cells= " << num_cells << std::endl;
  std::cout << "num_threads= " << num_threads << std::endl;

  // create a mt19937 instance for each thread
  for (int idx=0; idx<num_threads; idx++)
  {
    std::mt19937 gen(rd());
    mt19937_vec.push_back(gen);
  }

  #pragma omp parallel
  {
    #pragma omp single
    std::cout << "omp_get_num_threads() (parallel section) = " << omp_get_num_threads() << std::endl;
  }

  // omp_set_num_threads(8);
  // SeedRandom(0);

  // bool seed_per_thread = false;
  bool seed_per_thread = true;
  if (seed_per_thread)
    std::cout << "seed_per_thread = True" << std::endl;
  else
    std::cout << "seed_per_thread = False" << std::endl;

  if (seed_per_thread)
  {
    #pragma omp parallel private(num_threads, tid)
    {
      // int tid = omp_get_thread_num();
      tid = omp_get_thread_num();
      #pragma omp critical
      {
        // std::cout << "SeedRandom = " << tid << "\n" << std::flush;
        mt19937_vec[tid].seed(tid);
      }
    }
  }
  else
  {
    std::cout << "----- oops, same seed per thread is bad\n" << std::flush;
    exit(-1);
  }

  int chunksize = num_threads;
  // #pragma omp parallel for schedule(static,2)
  // #pragma omp parallel for schedule(static,4) private(rnum)
  // #pragma omp parallel for schedule(dynamic,chunksize) private(tid, rnum)
  // #pragma omp parallel for schedule(static,chunksize) private(tid, rnum)

  #pragma omp parallel for schedule(static) private(tid, rnum)
  for (int idx=0; idx<num_cells; idx++)
  {
    tid = omp_get_thread_num();
    rnum = UniformRandom(mt19937_vec[tid]);

    // #pragma omp critical  // otherwise, a " malloc: *** error for object 0x7fbe21600b48: pointer being freed was not allocated"
    // {
    //   rval.push_back(rnum);
    // }

    #pragma omp critical
    {
      // std::cout << "thread, idx, rnum = " << tid << ", "<< idx << ", " << rnum << "\n" << std::flush;
      // rval.push_back(rnum);   // Wrong
      rval[idx] = rnum;
    }
  }

  for (int idx=0; idx<num_cells; idx++)
      std::cout << "idx, rval = " << idx << ", " << rval[idx] << std::endl;
  return 0;
}