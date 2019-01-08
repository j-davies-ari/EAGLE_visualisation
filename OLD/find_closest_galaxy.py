import numpy as np
from sys import exit
import eagle as E

class groupnum_table(object):
    def __init__(self):
        ##############################################################################################################
        ref_sim = 'L0025N0376'
        ref_run = 'REFERENCE'
        sims = ['L0025N0376','L0025N0752']
        runs = [['StrongFB','WeakFB'],['REFERENCE','RECALIBRATED']]
        snap = '028_z000p000'
        infotype = 'SUBFIND'
        table = 'Subhalo'
        ##############################################################################################################

        simloc = '/data5/simulations/EAGLE/'+ref_sim+'/'+ref_run+'/data/'
        num_subs = np.array(E.readArray("SUBFIND_GROUP", simloc, snap, "/FOF/NumOfSubhalos"))
        masslist = np.array(E.readArray("SUBFIND_GROUP",simloc,snap,'FOF/Group_M_Crit200')[num_subs>0])*1e10
        COPlist = np.array(E.readArray("SUBFIND_GROUP",simloc,snap,'FOF/GroupCentreOfPotential')[num_subs>0])
        #gns = np.where(num_subs>0)[0]+1
        gns = np.arange(len(masslist)) + 1

        # Perform mass cut to remove dwarfs which pass close to COP
        masscut = np.where(np.log10(masslist)>11.7)[0]
        ref_COPlist = COPlist[masscut]
        self.ref_gns = gns[masscut]
        ref_masslist = masslist[masscut]

        self.connect_gns = {ref_sim+'_'+ref_run:self.ref_gns}

        for s, sim in enumerate(sims):
            for r, run in enumerate(runs[s]):

                simloc = '/data5/simulations/EAGLE/'+sim+'/'+run+'/data/'
                num_subs = np.array(E.readArray("SUBFIND_GROUP", simloc, snap, "/FOF/NumOfSubhalos"))
                masslist = np.array(E.readArray("SUBFIND_GROUP",simloc,snap,'FOF/Group_M_Crit200')[num_subs>0])*1e10
                COPlist = np.array(E.readArray("SUBFIND_GROUP",simloc,snap,'FOF/GroupCentreOfPotential')[num_subs>0])
                gns = np.where(num_subs>0)[0]+1

                masscut = np.where(np.log10(masslist)>11.7)[0]
                COPlist = COPlist[masscut]
                gns = gns[masscut]
                masslist = masslist[masscut]

                thisrun_equivalents = np.zeros(len(self.ref_gns))

                for n in range(len(self.ref_gns)):
                    COPdiff = np.sqrt((COPlist[:,0]-ref_COPlist[n,0])**2 + (COPlist[:,1]-ref_COPlist[n,1])**2 + (COPlist[:,2]-ref_COPlist[n,2])**2)
                    thisrun_equivalents[n] = gns[np.argmin(COPdiff)]
                self.connect_gns[sim+'_'+run] = thisrun_equivalents.astype(int)

    def matching_group(self,ref_group_number,other_sim,other_run):
        if ref_group_number not in self.ref_gns:
            raise ValueError('Group no. %s does not exist in Ref-L0025N0376 or is too low-mass (<10^11.7)'%(str(ref_group_number)))
        index = np.where(self.ref_gns==ref_group_number)[0]
        return self.connect_gns[other_sim+'_'+other_run][index][0]


if __name__ == '__main__':

    gn_table = groupnum_table()
    
    print gn_table.matching_group(32,'L0025N0376','StrongFB')
    print gn_table.matching_group(32,'L0025N0376','WeakFB')
    print gn_table.matching_group(32,'L0025N0752','REFERENCE')
    print gn_table.matching_group(32,'L0025N0752','RECALIBRATED')






