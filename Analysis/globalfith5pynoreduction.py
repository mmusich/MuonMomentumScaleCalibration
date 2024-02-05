import os
import multiprocessing
nthreads = multiprocessing.cpu_count()

os.environ["OMP_NUM_THREADS"] = f"{nthreads}"

import math
import ROOT
import numpy as np
import concurrent.futures
import matplotlib.pyplot as plt
import threading
import scipy
import scipy.linalg
import sys
import scipy.sparse.linalg
import scipy.sparse
import hdf5plugin
import h5py

np.set_printoptions(threshold=sys.maxsize)


def ptparm(parms):
  p = np.abs(1./parms[0])
  theta = np.pi/2. - parms[1]
  pt = p*np.sin(theta)
  return pt
  
def deltaphi(phi1,phi2):
  if phi1 < -np.pi:
    phi1 += 2*np.pi
  elif phi1 >= np.pi:
    phi1 += -2*np.pi

  if phi2 < -np.pi:
    phi2 += 2*np.pi
  elif phi2 >= np.pi:
    phi2 += -2*np.pi
  
  dphi = phi2-phi1
  if dphi < -np.pi:
    dphi += 2*np.pi
  elif dphi >= np.pi:
    dphi += -2*np.pi
  return dphi

filenamecor = "correctionResults_realigned.root"

filenameinfo = "/scratch/bruschin/cmsasymow/CMSSW_10_6_26/src/Analysis/HitAnalyzer/crab_submit/output/globalcor_data_0_1.root"
finfo = ROOT.TFile.Open(filenameinfo)

runtree = finfo.runtree

nparmsfull = np.int64(runtree.GetEntries())


oldcors = np.zeros((nparmsfull,), dtype=np.float64)

dooldcors = False

if dooldcors:
    fcor = ROOT.TFile.Open(filenamecor)
    parmtree = fcor.Get("parmtree")

    #oldcors = []
    for iparm,entry in enumerate(parmtree):
        oldcors[iparm] = parmtree.x
        #oldcors.append(parmtree.x)


oldcorstmp = np.zeros((nparmsfull,), dtype=np.float64)

dooldcorstmp = False

if dooldcorstmp:
    fcor = ROOT.TFile.Open(filenamecor)
    parmtree = fcor.Get("parmtree")

    #oldcors = []
    for iparm,entry in enumerate(parmtree):
        oldcorstmp[iparm] = parmtree.x
        #oldcors.append(parmtree.x)

testsize = nparmsfull


countsfull = np.zeros((testsize,), dtype=np.float64)
gradfull = np.zeros((testsize,1), dtype=np.float64)
hessfull = np.zeros((testsize, testsize), dtype=np.float64)

fgradnames.append("combinedgrads.hdf5")

weights = [1.]*len(fgradnames)

for weight, fgradname in zip(weights,fgradnames):
    print(fgradname, weight)
    with h5py.File(fgradname) as fgrads:
        gradd = fgrads["grad"]
        hessd = fgrads["hess"]
        # countsd = fgrads["counts"]
        
        # countsfull += countsd[...]
        
        gradfull[:,0] += weight*0.5*gradd[...]

        for i in range(nparmsfull):
            if i%5000 == 0:
                print(i)
                
            hessfull[i] += 0.5*np.abs(weight)*hessd[i]

parmset = set()
parmlistfull = []
for iidx,parm in enumerate(runtree):
    if iidx == testsize:
        break
    parmtype = runtree.parmtype
    ieta = math.floor(runtree.eta/0.1)
    iphi = math.floor(runtree.phi/(math.pi/8.))
    subdet = runtree.subdet
    layer = abs(runtree.layer)
    key = (parmtype, subdet, layer, (ieta,iphi))
    parmset.add(key)
    parmlistfull.append(key)

parmlist = parmlistfull

parmmap = {}
for iparm,key in enumerate(parmlist):
  parmmap[key] = iparm
  
nglobal = gradfull.shape[0]
print(nglobal)
grad = gradfull
hess = hessfull
idxmap = np.arange(nglobal)

isvalid = np.ones((nglobal,), dtype=np.bool_)

print("applying constraints and zeroing gradients for fixed parameters")

for iparm,parm in enumerate(parmlist):
    parmtype, subdet, layer, ieta = parm
    if parmtype==-1:
        print(f"null parameter for index {iparm}")
        grad[iparm] = 0.
        hess[iparm, :] = 0.
        hess[:, iparm] = 0.
        hess[iparm,iparm] = 2.
        
    siga = 0.
    if parmtype == 0:
        # local x translation
        siga = 5e-3

    elif parmtype == 1:
        # local y translation
        if subdet < 2:
            siga = 5e-3
        else:
            siga = 2.0

    elif parmtype == 2:
        siga = 5e-1
        #siga = 5e-3
    elif parmtype in [3, 4]:
        # out-of-plane rotations
        siga = 5e-3
    elif parmtype == 5:
        # in-plane rotation
        siga = 5e-3
    elif parmtype == 6:
        #b-field
        siga = 0.038
    elif parmtype == 7:
        # material energy loss
        siga = math.log(2.)

    elif parmtype in [8, 9]:
        #hit resolution (relative uncertainty on variance)
        #siga = 0.2
        siga = math.log(3.)

    elif parmtype == 10:
        # material resolution
        # siga = math.log(2.)
        siga = math.log(3.)
    elif parmtype == 11:
        # ionization resolution
        # siga = math.log(2.)
        siga = math.log(100.)

    #reco
    if False:
        grad[iparm] = 0.
        hess[iparm, :] = 0.
        hess[:,iparm] = 0.
        isvalid[iparm] = False

    oldhessdiag = hess[iparm, iparm]    
        
    if siga > 0.:
        grad[iparm] += oldcors[iparm]/siga**2
        hess[iparm, iparm] += 1./siga**2
    
    # sumsqoffdiag = np.sum(np.square(hess[iparm])) - np.square(hess[iparm,iparm])


    #if parmtype in [8] and oldhessdiag < 0.:
    if oldhessdiag < 0.:
      print("debug:", parmtype, subdet, layer, oldhessdiag, hess[iparm, iparm])
      sumsqoffdiag = np.sum(np.square(hess[iparm])) - np.square(hess[iparm,iparm])
      print("siga", siga)
      print("hessadd", 1./siga**2)
      print("sumsqoffdiag", sumsqoffdiag)


print("filling in lower triangular part of hessian")
def filllower(i):
  hess[i,:i] = hess[:i,i]
  
with concurrent.futures.ThreadPoolExecutor() as e:
  results = e.map(filllower, range(hess.shape[0]))
  
for result in results:
  pass

print("decomposing")

# transpose here is needed for the decomposition to occur in-place which speeds things up and saves a factor of 2 in memory
lower = True

chol = scipy.linalg.cho_factor(hess.transpose(), lower=lower, overwrite_a = True)

print("solving")
xout = scipy.linalg.cho_solve(chol, -grad)

print("done solve")

docov = False
doreplicas = False



if False:
    print("make identity")
    ident = np.identity(testsize, dtype = np.float64)
    print("compute inverse")
    cov = scipy.linalg.cho_solve(chol, ident.transpose(), overwrite_b = True)
    errs = np.sqrt(np.diag(cov))

    print("done compute inverse")

    neig = 500

    print("eigh")

    e, v = scipy.sparse.linalg.eigsh(cov, k=neig, which = "LA", tol = 0.1)

    # print("e", e)

    for ieig in range(neig):
        print("ieig", ieig)
        print(e[ieig])
        iv = v[:, ieig]
        imax = np.argmax(np.abs(iv))
        maxval = iv[imax]
        print("imax", imax)
        print("maxval", maxval)



    with h5py.File("reducedeig.hdf5", "w") as f:
        f.create_dataset("e", data = e, chunks=True, **hdf5plugin.Blosc(cname="lz4"))
        f.create_dataset("v", data = v, chunks=True, **hdf5plugin.Blosc(cname="lz4"))

    assert(0)

if docov:

    print("filling in upper triangular part of cholesky factor (with zeros)")
    def zero_upper_chol(i):
        chol[0][i,i+1:] = 0.

    with concurrent.futures.ThreadPoolExecutor() as e:
        results = e.map(zero_upper_chol, range(chol[0].shape[0]))


    print("computing cholinv")
    cholinv = scipy.linalg.lapack.dtrtri(chol[0], lower=chol[1], unitdiag = False, overwrite_c=True)[0]
    
    print("computing cov")
    cov = cholinv.T@cholinv
    
    print("computing errs")
    errs = np.sqrt(np.diag(cov))
    
    print("eigh")
    # e, v = scipy.sparse.linalg.eigsh(hess, k=1000, sigma = 0., OPinv = hess_solve_op)
    #
    e, v = scipy.sparse.linalg.eigsh(cov, k=10, which = "LM")

    print("e", e)

    assert(0)


    with h5py.File("correctionResults.hdf5", "w") as fout:
        print("writing h5py output")
        xouth5 = fout.create_dataset("x", (testsize,), dtype=np.float64, **hdf5plugin.LZ4())
        errsouth5 = fout.create_dataset("errs", (testsize,), dtype=np.float64, **hdf5plugin.LZ4())
        covouth5 = fout.create_dataset("cov", (testsize, testsize), dtype=np.float64, **hdf5plugin.LZ4(), chunks=(1, testsize))

        xouth5[...] = xout[:, 0]
        errsouth5[...] = errs
        for i in range(testsize):
            if i%1000 == 0:
                print(i)
            covouth5[i] = cov[i]

else:
    errs = np.zeros_like(np.diag(hess))





print("compute toys")

ntoys = 100



  
  
if doreplicas:
    print("computing toys")

    print("compute cov chol")

    cov_chol = scipy.linalg.cho_factor(cov, lower=True, overwrite_a = True)

    print("filling in upper triangular part of cholesky factor")
    def zero_upper(i):
      cov_chol[0][i,i+1:] = 0.

    with concurrent.futures.ThreadPoolExecutor(64) as e:
      results = e.map(zero_upper, range(cov_chol[0].shape[0]))

    for result in results:
      pass
  
    u = np.random.standard_normal((testsize,ntoys))

    print("xout.shape", xout.shape)
    xtoys = xout + cov_chol[0] @ u
    xtoys = np.where(isvalid[:, np.newaxis], xtoys, xout)
    print(xtoys.shape)



    print("done compute toys")



    errtoys = np.std(xtoys, axis = -1)

    err_ratio = errtoys/errs

    print(errs)
    print(errtoys)
    print(err_ratio)


#write output file
fout = ROOT.TFile.Open("correctionResults_lbl2018_recjpsidata.root", "RECREATE")

print("first loop")

idxmaptree = ROOT.TTree("idxmaptree","")
idx = np.empty((1), dtype=np.uint32)
idxmaptree.Branch("idx", idx, "idx/i")
for i in range(testsize):
    idx[0] = idxmap[i]
    idxmaptree.Fill()
    
idxmaptree.Write()

parmtree = ROOT.TTree("parmtree","")
x = np.empty((1), dtype=np.float32)
err = np.empty((1), dtype=np.float32)


xreplicas = np.empty((ntoys), dtype=np.float32)

print("second loop")

parmtree.Branch("x", x, "x/F")
parmtree.Branch("err", err, "err/F")
if doreplicas:
    parmtree.Branch("xreplicas", xreplicas , f"xreplicas[{ntoys}]/F")

for i in range(nglobal):
    x[0] = xout[i] + oldcors[i]
    err[0] = errs[i]
    if doreplicas:
        xreplicas[...] = xtoys[i]
    parmtree.Fill()
    
parmtree.Write()
fout.Close()
