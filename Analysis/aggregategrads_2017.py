import os

os.environ["XRD_PARALLELEVTLOOP"] = "24"


import ROOT

# ROOT.gInterpreter.ProcessLine(".O3")
#ROOT.ROOT.EnableImplicitMT(64)
ROOT.ROOT.EnableImplicitMT()


import numpy as np
import hdf5plugin
import h5py
# from utils import lumitools
import wremnants
import narf
import narf.lumitools


status = ROOT.gInterpreter.Declare('#include "aggregategrads.cpp"')

assert(status)

etas = np.linspace(-2.4, 2.4, 49)
logpts = np.linspace(-2.5,2.5,41)
cosphis = np.array([round(-1. + 2.*i/20,2) for i in range(21)],dtype=np.float64)
masses = np.linspace(2.9, 3.3, 101)


print(etas)
print(logpts)
print(cosphis)
print(masses)

neta = etas.shape[0] - 1
nlogpt = logpts.shape[0] - 1
ncosphi = cosphis.shape[0] - 1
nmass = masses.shape[0] - 1

fullidxs = -1*np.ones((neta+2,neta+2,nlogpt+2,ncosphi+2),dtype=np.int32)

datafits = False

if datafits:

    filenamefits = "fitsJDATA.hdf5"
    #ffits = h5py.file(filen
    with h5py.File(filenamefits) as ffits:
        print(ffits.keys())
        #assert(0)

        sigpdf = ffits["sigpdf"][...]
        bkgpdf = ffits["bkgpdf"][...]

        wsig = sigpdf/(sigpdf+bkgpdf)



        print("wsig.shape", wsig.shape)

        #assert(0)


        good_idx = ffits["good_idx"][...].astype(np.int32) + 1
        good_idx = (good_idx[0], good_idx[1], good_idx[2], good_idx[3])

        linidxs = np.arange(good_idx[0].shape[0])
        print(linidxs)

        #assert(0)

        #print(good_idx)

        fullidxs[good_idx] = linidxs



        #valid = np.where(fullidxs >= 0)

        #print(valid)

        #hdata = full[good_idx]

        #print(good_idx.shape)
        #print(fullidxs)
        #print(hdata)
        #print(full)


        #print(ffits["good_idx"])




    @ROOT.Numba.Declare(["float", "float", "float", "float", "float", "float", "float", "unsigned int"], "double")
    def massweight(ptplus, etaplus, phiplus, ptminus, etaminus, phiminus, mass, run):

        #if run == 1:
            #return 1.

        logpt = np.log(ptplus/ptminus)
        cosphi = np.cos(phiplus - phiminus)

        idx0 = np.digitize([etaplus], etas)
        idx1 = np.digitize([etaminus], etas)
        idx2 = np.digitize([logpt], logpts)
        idx3 = np.digitize([cosphi], cosphis)

        idxt = (idx0[0], idx1[0], idx2[0], idx3[0])

        #idxt = (np.digitize(etaplus, etas), np.digitize(etaminus, etas), np.digitize(logpt, logpts), np.digitize(cosphi, cosphis))

        idx = fullidxs[idxt]

        if idx == -1:
            return 0.

        massidx = np.digitize([mass], masses) - 1
        massidx = massidx[0]

        if massidx < 0 or massidx >= wsig.shape[-1]:
            return 0.


        if run == 1:
            return 1.

        return wsig[idx, massidx]

hlt_paths = ["HLT_DoubleMu4_3_Bs",
             "HLT_DoubleMu4_3_Jpsi_Displaced",
             "HLT_DoubleMu4_JpsiTrk_Displaced",
             "HLT_DoubleMu4_PsiPrimeTrk_Displaced",
             "HLT_Mu7p5_L2Mu2_Jpsi",
             "HLT_Mu7p5_Track2_Jpsi",
             "HLT_Mu7p5_Track3p5_Jpsi",
             "HLT_Mu7p5_Track7_Jpsi",
             "HLT_Dimuon0_LowMass_L1_0er1p5R",
             "HLT_Dimuon0_LowMass_L1_4R",
             "HLT_DoubleMu4_Jpsi_Displaced",
             "HLT_DoubleMu4_Jpsi_NoVertexing",
             "HLT_Dimuon10_PsiPrime_Barrel_Seagulls",
             "HLT_Dimuon20_Jpsi_Barrel_Seagulls",
             "HLT_Dimuon18_PsiPrime",
             "HLT_Dimuon25_Jpsi",
             "HLT_Dimuon0_Jpsi_L1_NoOS",
             "HLT_Dimuon0_Jpsi_NoVertexing_NoOS",
             "HLT_Dimuon0_Jpsi",
             "HLT_Dimuon0_Jpsi_NoVertexing",
             "HLT_Dimuon0_Jpsi_L1_4R_0er1p5R",
             "HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R",
             "HLT_Dimuon0_Jpsi3p5_Muon2",
             "HLT_Dimuon0_LowMass_L1_0er1p5",
             "HLT_Dimuon0_LowMass",
             "HLT_Dimuon0_LowMass_L1_4",
             "HLT_DoubleMu4_JpsiTrkTrk_Displaced",
             "HLT_Dimuon18_PsiPrime_noCorrL1",
             "HLT_Dimuon25_Jpsi_noCorrL1",
             "HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi",
             "HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi"]

chainjpsi = ROOT.TChain("tree")

chainjpsi.Add("/gpfs/ddn/srm/cms/store/user/musich/Charmonium/Charmonium2017B_v722_layerbylayer/240229_134706/0000/globalcor_data_0_*.root") #B

chainjpsi.Add("/gpfs/ddn/srm/cms/store/user/musich/Charmonium/Charmonium2017C_v722_layerbylayer/240229_134955/0000/globalcor_data_0_*.root") #C
chainjpsi.Add("/gpfs/ddn/srm/cms/store/user/musich/Charmonium/Charmonium2017C_v722_layerbylayer/240229_134955/0001/globalcor_data_0_*.root") #C
chainjpsi.Add("/gpfs/ddn/srm/cms/store/user/musich/Charmonium/Charmonium2017C_v722_layerbylayer/240229_134955/0002/globalcor_data_0_*.root") #C
chainjpsi.Add("/gpfs/ddn/srm/cms/store/user/musich/Charmonium/Charmonium2017C_v722_layerbylayer/240229_134955/0003/globalcor_data_0_*.root") #C
chainjpsi.Add("/gpfs/ddn/srm/cms/store/user/musich/Charmonium/Charmonium2017C_v722_layerbylayer/240229_134955/0004/globalcor_data_0_*.root") #C
chainjpsi.Add("/gpfs/ddn/srm/cms/store/user/musich/Charmonium/Charmonium2017C_v722_layerbylayer/240229_134955/0005/globalcor_data_0_*.root") #C

chainjpsi.Add("/gpfs/ddn/srm/cms/store/user/musich/Charmonium/Charmonium2017D_v722_layerbylayer/240229_135157/0000/globalcor_data_0_*.root") #D
chainjpsi.Add("/gpfs/ddn/srm/cms/store/user/musich/Charmonium/Charmonium2017D_v722_layerbylayer/240229_135157/0001/globalcor_data_0_*.root") #D
chainjpsi.Add("/gpfs/ddn/srm/cms/store/user/musich/Charmonium/Charmonium2017D_v722_layerbylayer/240229_135157/0002/globalcor_data_0_*.root") #D

chainjpsi.Add("/gpfs/ddn/srm/cms/store/user/musich/Charmonium/Charmonium2017E_v722_layerbylayer/240229_135435/0000/globalcor_data_0_*.root") #E
chainjpsi.Add("/gpfs/ddn/srm/cms/store/user/musich/Charmonium/Charmonium2017E_v722_layerbylayer/240229_135435/0001/globalcor_data_0_*.root") #E
chainjpsi.Add("/gpfs/ddn/srm/cms/store/user/musich/Charmonium/Charmonium2017E_v722_layerbylayer/240229_135435/0002/globalcor_data_0_*.root") #E
chainjpsi.Add("/gpfs/ddn/srm/cms/store/user/musich/Charmonium/Charmonium2017E_v722_layerbylayer/240229_135435/0003/globalcor_data_0_*.root") #E

chainjpsi.Add("/gpfs/ddn/srm/cms/store/user/musich/Charmonium/Charmonium2017F_v722_layerbylayer/240229_135605/0000/globalcor_data_0_*.root") #F
chainjpsi.Add("/gpfs/ddn/srm/cms/store/user/musich/Charmonium/Charmonium2017F_v722_layerbylayer/240229_135605/0001/globalcor_data_0_*.root") #F
chainjpsi.Add("/gpfs/ddn/srm/cms/store/user/musich/Charmonium/Charmonium2017F_v722_layerbylayer/240229_135605/0002/globalcor_data_0_*.root") #F
chainjpsi.Add("/gpfs/ddn/srm/cms/store/user/musich/Charmonium/Charmonium2017F_v722_layerbylayer/240229_135605/0003/globalcor_data_0_*.root") #F
chainjpsi.Add("/gpfs/ddn/srm/cms/store/user/musich/Charmonium/Charmonium2017F_v722_layerbylayer/240229_135605/0004/globalcor_data_0_*.root") #F

wremdir = os.environ["WREM_BASE"]

jsonhelper = narf.lumitools.make_jsonhelper(f"{wremdir}/wremnants-data/data/Cert_294927-306462_13TeV_UL2017_Collisions17_HLT_IsoMu24_v_CustomJSON.txt")

filenameinfo = chainjpsi.GetListOfFiles()[0].GetTitle()
finfo = ROOT.TFile.Open(filenameinfo)
runtree = finfo.Get("runtree")
nparms = int(runtree.GetEntries())

dj = ROOT.ROOT.RDataFrame(chainjpsi)


print(dj.Sum("edmvalref").GetValue())
# print(dj.Max("gradmax").GetValue())

# assert(0)

dj = dj.Filter(jsonhelper, ["run", "lumi"], "jsonfilter")

dj = dj.Filter(" || ".join(hlt_paths))

dj = dj.Filter("Muplus_nvalid > 5 && Muplus_nvalidpixel>0 && Muminus_nvalid > 5 && Muminus_nvalidpixel > 0 && Jpsi_mass > 2.9 && Jpsi_mass < 3.3")

dj = dj.Filter("Mupluscons_pt > 5.0 && Muminuscons_pt > 5.0")

dj = dj.Filter("edmvalref < 1e-5")

dj = dj.Filter("chisqval/ndof < 3.")

#dj = dj.Define("massweightval", "Numba::massweight(Muplus_pt, Muplus_eta, Muplus_phi, Muminus_pt, Muminus_eta, Muminus_phi, Jpsi_mass, run)")

dj = dj.Define("massweightval", "1.0")

dj = dj.Filter("massweightval > 0.")


djcount = dj.Count()
maxgradient = dj.Max("gradmax")


debugweights = False

if debugweights:
    hmass = dj.Histo1D(("hmass", "", 100, 2.9, 3.3), "Jpsi_mass")
    hmassweighted = dj.Histo1D(("hmassweighted", "", 100, 2.9, 3.3), "Jpsi_mass", "massweightval")


    c = ROOT.TCanvas()
    hmass.SetLineColor(ROOT.kRed)
    hmass.Draw("HIST")
    hmassweighted.Draw("HISTSAME")


    input("weight")

gradhelperj = ROOT.GradHelper(nparms)
hesshelperj = ROOT.HessHelper(nparms)

grad = dj.Book(gradhelperj, ["gradv", "globalidxv", "massweightval"])
hess = dj.Book(hesshelperj, ["hesspackedv", "globalidxv", "massweightval"])


#gradval = grad.GetResult()
#hessval = hess.GetResult()

#print(gradval[0])

print("djcount", djcount.GetValue())
print("maxgradient", maxgradient.GetValue())


#chunksize = 32

#fout = h5py.File("combinedgrads.hdf5", "w", rdcc_nbytes = nparms*8*chunksize*4, rdcc_nslots = nparms//chunksize*10)

fout = h5py.File("combinedgrads.hdf5", "w")

#gradout = fout.create_dataset("grad", (nparms,), dtype=np.float64, compression="lzf")
#hessout = fout.create_dataset("hess", (nparms, nparms), dtype=np.float64, compression="lzf", chunks=(1, nparms))

gradout = fout.create_dataset("grad", (nparms,), dtype=np.float64, **hdf5plugin.LZ4())
hessout = fout.create_dataset("hess", (nparms, nparms), dtype=np.float64, **hdf5plugin.LZ4(), chunks=(1, nparms))



gradout[...] = grad

#assert(0)

hessrow = np.zeros((nparms,), dtype=np.float64)

for i in range(nparms):
  if i%1000 == 0:
      print(i)
  #if i>61200:
    #print(i)
  hess.fill_row(i, hessrow)
  hessout[i] = hessrow
  #hess.fill_row(i, hessout[i])
  #if i < 10:
    #print(hessout[i])

#hess.fill_row(0, hessrow.data)

#print(hessrow)
