from pymongo import MongoClient
import sys
import petaloAnalysis

#connect to db
client = MongoClient('127.0.0.1', 27017)
db = client['petalo']
collection = db['cherenkov']

#load data
petalo = petaloAnalysis.PetaloAnalysis('sensors.csv','events.csv')

def processData(params):
    vel = params[0]
    pde = params[1]
    minwl = params[2]
    maxwl = params[3]
    asic = params[4]
    sptr = params[5]
    npe = params[6]
    print "vel: {}, pde: {}, asic: {}, sptr: {}, minwl: {}, maxwl: {}, npe: {}".format(vel,pde,asic,sptr,minwl,maxwl,npe)

    petalo.processWaveform(minwl,maxwl,pde,asic,sptr)
    result = petalo.computeAVG(vel,npe)

    data = {}
    data['vel'] = vel
    data['pde'] = pde
    data['asic'] = asic
    data['sptr'] = sptr
    data['minwl'] = minwl
    data['maxwl'] = maxwl
    data['npe'] = npe

    data['avg'] = result.tolist()

    collection.insert_one(data)



if __name__ == "__main__":
    coreNumber = int(sys.argv[1])
    totalCores = int(sys.argv[2])

    pdes = [x for x in xrange(10,105,5)]
    minwls = [x for x in xrange(150,500,50)]
    asics = [x for x in xrange(0,40,10)]
    sptrs = [x for x in xrange(0,90,10)]

    configs = [(210,pde,minwl,1500, asic, sptr,10) for pde in pdes for minwl in minwls for asic in asics for sptr in sptrs]

    jobsPerCore = len(configs)/totalCores
    start = coreNumber*jobsPerCore
    if coreNumber < totalCores-1:
        end = start + jobsPerCore
    else:
        end = len(configs)

    configsThisCore = configs[start:end]

    map(processData, configsThisCore)
