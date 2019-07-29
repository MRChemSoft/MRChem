#pragma once

#include "MRCPP/Parallel"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include "mrchem.h"
#include "qmfunctions/qmfunction_fwd.h"
#include "utils/Bank.h"

namespace mrchem {

namespace omp {
extern int n_threads;
} // namespace omp

namespace mpi {
extern bool numerically_exact;
extern int shared_memory_size;

extern int world_rank;
extern int world_size;
extern int orb_rank;
extern int orb_size;
extern int share_rank;
extern int share_size;
extern int sh_group_rank;
extern int is_bank;
extern int is_bankclient;
extern int bankmaster;

extern MPI_Comm comm_orb;
extern MPI_Comm comm_share;
extern MPI_Comm comm_sh_group;
extern MPI_Comm comm_bank;

extern Bank orb_bank;

void initialize();
void finalize();
void barrier(MPI_Comm comm);

bool grand_master();
bool share_master();
bool my_orb(const Orbital &orb);
bool my_unique_orb(const Orbital &orb);
void distribute(OrbitalVector &Phi);
void free_foreign(OrbitalVector &Phi);
OrbitalChunk get_my_chunk(OrbitalVector &Phi);

void send_orbital(Orbital &orb, int dst, int tag, MPI_Comm comm = mpi::comm_orb);
void recv_orbital(Orbital &orb, int src, int tag);

void send_function(QMFunction &func, int dst, int tag, MPI_Comm comm);
void recv_function(QMFunction &func, int src, int tag, MPI_Comm comm);
void share_function(QMFunction &func, int src, int tag, MPI_Comm comm);

void reduce_function(double prec, QMFunction &func, MPI_Comm comm);
void broadcast_function(QMFunction &func, MPI_Comm comm);

void allreduce_vector(IntVector &vec, MPI_Comm comm);
void allreduce_vector(DoubleVector &vec, MPI_Comm comm);
void allreduce_vector(ComplexVector &vec, MPI_Comm comm);
void allreduce_matrix(IntMatrix &vec, MPI_Comm comm);
void allreduce_matrix(DoubleMatrix &mat, MPI_Comm comm);
void allreduce_matrix(ComplexMatrix &mat, MPI_Comm comm);

} // namespace mpi

} // namespace mrchem
