// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// This class is heavily inspired from the iMultiFab class in AmReX. Notable
// additions include the geometric shifting of the interface across periodic
// boundaries during parallel copy

#ifndef IRL_SEPUMULTIFAB_TPP_
#define IRL_SEPUMULTIFAB_TPP_

namespace amrex {

namespace detail {
void read_fab(SepUnionArrayBox& fab, VisMF::FabOnDisk const& fod,
              std::string const& name) {
  std::string fullname = VisMF::DirName(name);
  fullname += fod.m_name;
  VisMFBuffer::IO_Buffer io_buffer(VisMFBuffer::GetIOBufferSize());
  std::ifstream ifs;
  ifs.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
  ifs.open(fullname.c_str(), std::ios::in | std::ios::binary);
  if (!ifs.good()) {
    amrex::FileOpenFailed(fullname);
  }
  ifs.seekg(fod.m_head, std::ios::beg);
  fab.readFrom(ifs);
}
}  // namespace detail

void SepUnionMultiFab::ParallelCopyWithPeriodicShift(
    const SepUnionMultiFab& src, int scomp, int dcomp, int ncomp,
    const IntVect& snghost, const IntVect& dnghost, const Geometry& geom,
    CpOp op, const FabArrayBase::CPC* a_cpc, bool deterministic) {
  this->ParallelCopy(src, scomp, dcomp, ncomp, snghost, dnghost,
                     geom.periodicity(), op, a_cpc, deterministic);
  this->PeriodicShift(geom);
}

void SepUnionMultiFab::ParallelCopyWithPeriodicShift(
    const SepUnionMultiFab& src, int src_comp, int dest_comp, int num_comp,
    int src_nghost, int dst_nghost, const Geometry& geom) {
  this->ParallelCopy(src, src_comp, dest_comp, num_comp, src_nghost, dst_nghost,
                     geom.periodicity());
  this->PeriodicShift(geom);
}

void SepUnionMultiFab::FillBoundaryWithPeriodicShift(int scomp, int ncomp,
                                                     const Geometry& geom,
                                                     bool cross) {
  this->FillBoundary(scomp, ncomp, geom.periodicity(), cross);
  this->PeriodicShift(geom);
}

void SepUnionMultiFab::FillBoundaryWithPeriodicShift(const Geometry& geom) {
  this->FillBoundary(geom.periodicity());
  this->PeriodicShift(geom);
}

void SepUnionMultiFab::PeriodicShift(const Geometry& geom) {
  const Box& domain = geom.Domain();
  for (MFIter mfi((*this), TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    const Box& grown_box = mfi.growntilebox();           // valid + all ghosts
    const Box domain_intersection = grown_box & domain;  // cells inside domain
    const Array4<IRL::SeparatorUnion>& arr = (*this).array(mfi);

    // Loop over all ghost cells (grown_box \ domain_intersection)
    amrex::ParallelFor(grown_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      // Skip if inside the domain (not a ghost)
      if (domain_intersection.contains(i, j, k)) return;
      // Is my ghost periodic?
      bool is_periodic_ghost_x0 = false;
      bool is_periodic_ghost_x1 = false;
      bool is_periodic_ghost_y0 = false;
      bool is_periodic_ghost_y1 = false;
      bool is_periodic_ghost_z0 = false;
      bool is_periodic_ghost_z1 = false;
      AMREX_D_TERM(
          is_periodic_ghost_x0 = (geom.isPeriodic(0) && i < domain.smallEnd(0));
          is_periodic_ghost_x1 = (geom.isPeriodic(0) && i > domain.bigEnd(0));
          ,
          is_periodic_ghost_y0 = (geom.isPeriodic(1) && j < domain.smallEnd(1));
          is_periodic_ghost_y1 = (geom.isPeriodic(1) && j > domain.bigEnd(1));
          ,
          is_periodic_ghost_z0 = (geom.isPeriodic(2) && k < domain.smallEnd(2));
          is_periodic_ghost_z1 = (geom.isPeriodic(2) && k > domain.bigEnd(2)));
      const bool is_periodic_ghost =
          is_periodic_ghost_x0 || is_periodic_ghost_y0 ||
          is_periodic_ghost_z0 || is_periodic_ghost_x1 ||
          is_periodic_ghost_y1 || is_periodic_ghost_z1;
      if (!is_periodic_ghost) return;
      // Shift interface datum
      const IRL::Pt shift =
          IRL::Pt(is_periodic_ghost_x0
                      ? -geom.ProbLength(0)
                      : (is_periodic_ghost_x1 ? geom.ProbLength(0) : 0.0),
                  is_periodic_ghost_y0
                      ? -geom.ProbLength(1)
                      : (is_periodic_ghost_y1 ? geom.ProbLength(1) : 0.0),
                  is_periodic_ghost_z0
                      ? -geom.ProbLength(2)
                      : (is_periodic_ghost_z1 ? geom.ProbLength(2) : 0.0));
      arr(i, j, k).shift(shift);
    });
  }
}

void SepUnionMultiFab_Write(const SepUnionMultiFab& mf, const std::string& name,
                            bool set_ghost) {
  BL_PROFILE("Write(SepUnionMultiFab)");
  AMREX_ASSERT(name.back() != '/');

  int data_bytes = sizeof(IRL::SeparatorUnion);

  bool useSparseFPP = false;
  const Vector<int>& pmap = mf.DistributionMap().ProcessorMap();
  std::set<int> procsWithData;
  Vector<int> procsWithDataVector;
  for (int i : pmap) {
    procsWithData.insert(i);
  }
  const int nOutFiles = VisMF::GetNOutFiles();
  if (static_cast<int>(procsWithData.size()) < nOutFiles) {
    useSparseFPP = true;
    for (auto i : procsWithData) {
      procsWithDataVector.push_back(i);
    }
  }

  std::string filePrefix = name + "_D_";

  NFilesIter nfi(nOutFiles, filePrefix, VisMF::GetGroupSets(),
                 VisMF::GetSetBuf());

  if (useSparseFPP) {
    nfi.SetSparseFPP(procsWithDataVector);
  } else {
    nfi.SetDynamic();
  }

  const auto& fio = SepUnionArrayBox::getFABio();

  for (; nfi.ReadyToWrite(); ++nfi) {
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
      const SepUnionArrayBox& fab = mf[mfi];
      {
        std::stringstream hss;
        fio.write_header(hss, fab, fab.nComp());
        auto hLength = static_cast<std::streamoff>(hss.tellp());
        auto tstr = hss.str();
        nfi.Stream().write(tstr.c_str(), hLength);
        nfi.Stream().flush();
      }
      auto const* fabdata = fab.dataPtr();
#ifdef AMREX_USE_GPU
      std::unique_ptr<SepUnionArrayBox> hostfab;
      if (fab.arena()->isManaged() || fab.arena()->isDevice()) {
        hostfab = std::make_unique<SepUnionArrayBox>(fab.box(), fab.nComp(),
                                                     The_Pinned_Arena());
        Gpu::dtoh_memcpy_async(
            hostfab->dataPtr(), fab.dataPtr(),
            fab.size() * sizeof(typename SepUnionArrayBox::value_type));
        Gpu::streamSynchronize();
        fabdata = hostfab->dataPtr();
      }
#endif
      Long writeDataItems = fab.box().numPts() * mf.nComp();
      Long writeDataSize = writeDataItems * data_bytes;
      nfi.Stream().write((char*)fabdata, writeDataSize);
      nfi.Stream().flush();
    }
  }

  int coordinatorProc = ParallelDescriptor::IOProcessorNumber();
  if (nfi.GetDynamic()) {
    coordinatorProc = nfi.CoordinatorProc();
  }

  if (ParallelDescriptor::MyProc() == coordinatorProc) {
    std::string header_file_name = name + "_H";
    VisMFBuffer::IO_Buffer io_buffer(VisMFBuffer::GetIOBufferSize());
    std::ofstream ofs;
    ofs.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
    ofs.open(header_file_name.c_str(), std::ios::out | std::ios::trunc);
    if (!ofs.good()) {
      amrex::FileOpenFailed(header_file_name);
    }

    ofs << "amrex::FabArray<" << SepUnionArrayBox::getClassName() << "> v1.0\n";
    ofs << mf.nComp() << '\n';
    ofs << mf.nGrowVect() << '\n';
    mf.boxArray().writeOn(ofs);
    ofs << '\n';

    const DistributionMapping& dm = mf.DistributionMap();
    int nfabs = mf.boxArray().size();
    int nFiles = NFilesIter::ActualNFiles(nOutFiles);
    int nprocs = ParallelDescriptor::NProcs();

    Vector<Long> fabBytes(nfabs, 0);
    std::map<int, Vector<int>> rankBoxOrder;
    for (int i = 0; i < nfabs; ++i) {
      std::stringstream hss;
      SepUnionArrayBox tmp(mf.fabbox(i), mf.nComp(), false);
      fio.write_header(hss, tmp, tmp.nComp());
      // Size includes header and data
      fabBytes[i] =
          static_cast<std::streamoff>(hss.tellp()) + tmp.size() * data_bytes;
      rankBoxOrder[dm[i]].push_back(i);
    }

    Vector<int> fileNumbers;
    if (nfi.GetDynamic()) {
      fileNumbers = nfi.FileNumbersWritten();
    } else if (nfi.GetSparseFPP()) {
      // if sparse, write to (file number = rank)
      fileNumbers.resize(nprocs);
      std::iota(fileNumbers.begin(), fileNumbers.end(), 0);
    } else {
      fileNumbers.resize(nprocs);
      for (int i = 0; i < nprocs; ++i) {
        fileNumbers[i] =
            NFilesIter::FileNumber(nFiles, i, VisMF::GetGroupSets());
      }
    }

    Vector<VisMF::FabOnDisk> fod(nfabs);

    const Vector<Vector<int>>& fileNumbersWriteOrder =
        nfi.FileNumbersWriteOrder();
    for (auto const& rv : fileNumbersWriteOrder) {
      Long currentOffset = 0;
      for (auto rank : rv) {
        auto rbo_it = rankBoxOrder.find(rank);
        if (rbo_it != rankBoxOrder.end()) {
          Vector<int> const& index = rbo_it->second;
          int whichFileNumber = fileNumbers[rank];
          std::string const& whichFileName = VisMF::BaseName(
              NFilesIter::FileName(whichFileNumber, filePrefix));
          for (int i : index) {
            fod[i].m_name = whichFileName;
            fod[i].m_head = currentOffset;
            currentOffset += fabBytes[i];
          }
        }
      }
    }
    ofs << fod;
  }
}

void SepUnionMultiFab_Read(SepUnionMultiFab& mf, const std::string& name,
                           const char* faHeader, int coordinatorProc,
                           int allow_empty_mf) {
  BL_PROFILE("Read(SepUnionMultiFab)");
  AMREX_ASSERT(name.back() != '/');

  BoxArray ba;
  int ncomp;
  IntVect ngrow;
  Vector<VisMF::FabOnDisk> fod;
  {
    std::string header_file_name = name + "_H";
    Vector<char> header_file_chars;
    ParallelDescriptor::ReadAndBcastFile(header_file_name, header_file_chars);
    std::string header_file_string(header_file_chars.data());
    std::stringstream ifs(header_file_string, std::istringstream::in);

    std::string type, version;
    ifs >> type >> version;
    AMREX_ASSERT(type == "amrex::FabArray<amrex::SepUnionArrayBox>" ||
                 type == "amrex::FabArray<amrex::SepUnionArrayBox>");
    ifs >> ncomp;
    ifs >> ngrow;
    ba.readFrom(ifs);
    ifs >> fod;
  }

  if (mf.empty()) {
    mf.define(ba, DistributionMapping{ba}, ncomp, ngrow);
  } else {
    AMREX_ASSERT(amrex::match(ba, mf.boxArray()));
  }

#ifdef AMREX_USE_MPI
  const int nopensperfile =
      VisMF::GetMFFileInStreams();  // # of concurrent readers per file
  const int myproc = ParallelDescriptor::MyProc();
  const int coordproc = ParallelDescriptor::IOProcessorNumber();

  int nreqs = 0;
  int allreadsindex = 0;
  std::map<std::string, int> filenames;  // <filename, allreadsindex>

  const int nboxes = mf.size();
  const auto& dm = mf.DistributionMap();
  for (int i = 0; i < nboxes; ++i) {
    if (myproc == dm[i]) {
      ++nreqs;
    }
    if (myproc == coordproc) {
      std::string const& fname = fod[i].m_name;
      auto r = filenames.insert(std::make_pair(fname, allreadsindex));
      if (r.second) {
        ++allreadsindex;
      }
    }
  }

  const int readtag = ParallelDescriptor::SeqNum();
  const int donetag = ParallelDescriptor::SeqNum();

  if (myproc == coordproc) {
    std::multiset<int>
        availablefiles;  // [whichFile]  supports multiple reads/file
    Vector<std::map<int, std::map<Long, int>>>
        allreads;  // [file]<proc,<seek,index>>

    const auto nfiles = static_cast<int>(filenames.size());
    for (int i = 0; i < nfiles; ++i) {
      for (int j = 0; j < nopensperfile; ++j) {
        availablefiles.insert(i);
      }
    }
    allreads.resize(nfiles);
    for (int i = 0; i < nboxes; ++i) {
      const auto whichproc = dm[i];
      const auto iseekpos = fod[i].m_head;
      std::string const& fname = fod[i].m_name;
      auto filenamesiter = filenames.find(fname);
      if (filenamesiter != filenames.end()) {
        const int fi = filenamesiter->second;
        allreads[fi][whichproc].insert(std::make_pair(iseekpos, i));
      } else {
        amrex::Error("Error in amrex::Read: filename not found " + fname);
      }
    }

    int totalioreqs = nboxes;
    int reqspending = 0;
    int iopfileindex;
    std::deque<int> iopreads;
    std::set<int> busyprocs;
    while (totalioreqs > 0) {
      auto afilesiter = availablefiles.begin();
      while (afilesiter != availablefiles.end()) {
        const int fi = *afilesiter;
        if (allreads[fi].empty()) {
          availablefiles.erase(fi);
          afilesiter = availablefiles.begin();
          continue;
        }
        auto whichread = allreads[fi].begin();
        for (; whichread != allreads[fi].end(); ++whichread) {
          const int tryproc = whichread->first;
          if (busyprocs.find(tryproc) == busyprocs.end()) {  // not busy
            busyprocs.insert(tryproc);
            Vector<int> vreads;
            vreads.reserve(whichread->second.size());
            for (auto const& kv : whichread->second) {
              vreads.push_back(kv.second);
            }
            if (tryproc == coordproc) {
              iopfileindex = fi;
              for (auto x : vreads) {
                iopreads.push_back(x);
              }
            } else {
              ParallelDescriptor::Send(vreads, tryproc, readtag);
              ++reqspending;
            }
            availablefiles.erase(afilesiter);
            afilesiter = availablefiles.begin();
            break;
          }
        }
        if (whichread == allreads[fi].end()) {
          ++afilesiter;
        } else {
          allreads[fi].erase(whichread);
        }
      }

      while (!iopreads.empty()) {
        int i = iopreads.front();
        detail::read_fab(mf[i], fod[i], name);
        --totalioreqs;
        iopreads.pop_front();
        if (iopreads.empty()) {
          availablefiles.insert(iopfileindex);
          busyprocs.erase(coordproc);
        }
        int doneflag;
        MPI_Status status;
        ParallelDescriptor::IProbe(MPI_ANY_SOURCE, donetag, doneflag, status);
        if (doneflag) {
          break;
        }
      }

      if (reqspending > 0) {
        Vector<int> idone(2);
        ParallelDescriptor::Message rmess =
            ParallelDescriptor::Recv(idone, MPI_ANY_SOURCE, donetag);
        const int i = idone[0];
        const int procdone = rmess.pid();
        totalioreqs -= idone[1];
        --reqspending;
        busyprocs.erase(procdone);
        std::string const& fname = fod[i].m_name;
        const int fi = filenames.find(fname)->second;
        availablefiles.insert(fi);
      }
    }
  } else {
    Vector<int> recreads(nreqs, -1);
    Vector<int> idone(2);
    while (nreqs > 0) {
      ParallelDescriptor::Message rmess =
          ParallelDescriptor::Recv(recreads, coordproc, readtag);
      const auto nrmess = static_cast<int>(rmess.count());
      for (int ir = 0; ir < nrmess; ++ir) {
        int i = recreads[ir];
        detail::read_fab(mf[i], fod[i], name);
      }
      nreqs -= nrmess;
      idone[0] = recreads[0];
      idone[1] = nrmess;
      ParallelDescriptor::Send(idone, coordproc, donetag);
    }
  }
#else
  for (MFIter mfi(fa); mfi.isValid(); ++mfi) {
    detail::read_fab(fa[mfi], fod[mfi.index()], name);
  }
#endif
}

}  // namespace amrex

#endif  // IRL_SEPUMULTIFAB_TPP_
