## This is a object-based version o Man Liang & Vaughan Voller's discritized
## delta xsectional model, after man_delta4.py.
## DEJ Hobley, Aug 2015.

# assumptions:
# The topset always maintains its slope - eroding if required
# In a rising BL the foreset can be abandoned
# Tank has constant level basement
# erosion w/o depo is always forbidden in a single section

# Most likely cause of this code crashing is it reaching the RHS of its domain

import numpy as np
from pylab import figure, plot, show, close
import matplotlib.pyplot as plt

class delta(object):
    """
    """

    def execute(self, input_file, SL_trajectory,
                print_strata=True, graphs=True, completeness_records=[],
                evolving_pys = False,
                restricted_channel_mass_conserved = False,
                decoupled_erosion_depo = False,
                walking_erosion_depo=False,
                compensation=False,
                never_erosion=False):
        # load the input params from the textfile:
        input_dict = self.read_input_file(input_file)
        n = int(input_dict['n'])
        delr = float(input_dict['delr'])  # note this is overwritten below
        delt = float(input_dict['delt'])
        nt = int(input_dict['nt'])
        self.delt = delt
        self.nt = nt
        Q = float(input_dict['Q'])
        ST = float(input_dict['ST'])
        SF = float(input_dict['SF'])
        theta = float(input_dict['theta'])

        rnode = np.arange(0.5, n+0.5, 1.)**0.5 * delr
        px_left_edge = np.zeros_like(rnode)
        px_left_edge[1:] = (rnode[1:] - rnode[:-1])/2. + rnode[1:]
        delr = rnode[1:]-rnode[:-1] #this is cludgey, should do accurately
        rnode = rnode[:-1]
        etanew = np.zeros_like(rnode)  # sed height in column
        eta = np.zeros_like(rnode)
        sedflux_modifier = 1.
        erosion_modifier = 1.

        if compensation:
            erosion_py_width = float(input_dict['erosion_py_width'])
            depo_py_width = float(input_dict['depo_py_width'])
            assert erosion_py_width <= depo_py_width
            sedflux_modifier /= depo_py_width
            erosion_modifier = erosion_py_width/depo_py_width  # as above
            section_position = np.random.rand()
            # ^these need to get halved for comparison with pys above
            drift = float(input_dict['drift'])
            pos_or_neg = 1.  #make this -1. to reverse the drift direction
            channel_position = 2.*(np.random.rand(nt)-0.5)*drift  # the steps
            channel_start = np.random.rand()
            channel_position[0] = channel_start
            depo_in_section = np.fabs(channel_start-section_position)<(
                                                           erosion_py_width/2.)
            # ^ this tells us if we're counting the init pos as "active" or
            # "inactive" depo, and will govern how and when we jump the depo
            # locus.
            # \/ this stuff needs to go into the active loop now
            for i in xrange(nt-1):
                channel_position[i+1] += pos_or_neg*channel_position[i]
                excess = channel_position[i+1]%1.
                if channel_position[i+1] > 1.:  # reverse direction from 1
                    pos_or_neg *= -1.
                    #print channel_position[i+1]
                    channel_position[i+1] = 1. - excess
                    #print channel_position[i+1]
                    #print '***'
                elif channel_position[i+1] < 0.:  # reverse direction from 0.
                    pos_or_neg *= -1.
                    channel_position[i+1] = 1. - excess
            historic_high = 0.
        # now the stuff conditional on the param rules:
        if evolving_pys:
            assert compensation is False, "Cannot use evolving_pys with compensation!!"
            if not rule_override:
                rule = int(input_dict['rule'])
            if rule == 1:
                using_rule_1 = True  # have to do this dynamically in the loop
                using_rule_2 = False
                historic_high = 0.
                channelised = False
            elif rule == 2:
                using_rule_1 = False
                using_rule_2 = True
            else:
                raise ValueError
        else:
            using_rule_1 = False
            using_rule_2 = False
        if restricted_channel_mass_conserved:
            activity_py = float(input_dict['activity_py'])
            erosion_py = activity_py
            depo_py = activity_py
            sedflux_modifier /= depo_py
        if decoupled_erosion_depo:
            erosion_py = float(input_dict['erosion_py'])
            depo_py = float(input_dict['depo_py'])
            assert erosion_py <= depo_py
            # ^we need to assume this to handle mass balance
            sedflux_modifier /= depo_py
            erosion_modifier = erosion_py/depo_py
        if walking_erosion_depo:
            erosion_py_width = float(input_dict['erosion_py_width'])
            depo_py_width = float(input_dict['depo_py_width'])
            assert erosion_py_width <= depo_py_width
            sedflux_modifier /= depo_py_width
            erosion_modifier = erosion_py_width/depo_py_width  # as above
            section_position = np.random.rand()
            # ^these need to get halved for comparison with pys above
            drift = float(input_dict['drift'])
            pos_or_neg = 1.  #make this -1. to reverse the drift direction
            channel_position = 2.*(np.random.rand(nt)-0.5)*drift  # the steps
            channel_start = np.random.rand()
            channel_position[0] = channel_start
            for i in xrange(nt-1):
                channel_position[i+1] += pos_or_neg*channel_position[i]
                excess = channel_position[i+1]%1.
                if channel_position[i+1] > 1.:  # reverse direction from 1
                    pos_or_neg *= -1.
                    #print channel_position[i+1]
                    channel_position[i+1] = 1. - excess
                    #print channel_position[i+1]
                    #print '***'
                elif channel_position[i+1] < 0.:  # reverse direction from 0.
                    pos_or_neg *= -1.
                    channel_position[i+1] = 1. - excess
        # all these methods need this...
        py_thresholds_depo = np.random.rand(nt)
        py_thresholds_inc = np.random.rand(nt)
        eroded_last_time = False

        #handling to speed up processing in the loop:
        modified_rules = np.any((restricted_channel_mass_conserved, never_erosion,
                                decoupled_erosion_depo, walking_erosion_depo,
                                compensation))
        assert np.sum((restricted_channel_mass_conserved, never_erosion,
                      decoupled_erosion_depo, walking_erosion_depo,
                      compensation)) <= 1
        if modified_rules and not evolving_pys:
            if restricted_channel_mass_conserved:
                stepwise_depo = np.greater(py_thresholds_depo, 1.-depo_py)
                stepwise_inc = stepwise_depo
            elif decoupled_erosion_depo:
                stepwise_depo = np.greater(py_thresholds_depo, 1.-depo_py)
                stepwise_inc = np.zeros_like(stepwise_depo)
                stepwise_inc[stepwise_depo] = np.greater(py_thresholds_inc[stepwise_depo], 1.-erosion_py/depo_py)
            elif never_erosion:
                stepwise_depo = np.ones(nt, dtype=bool)
                stepwise_inc = np.zeros(nt, dtype=bool)
            elif walking_erosion_depo or compensation:
                stepwise_depo = np.less(np.absolute(section_position-channel_position),
                                        depo_py_width)
                stepwise_inc = np.less(np.absolute(section_position-channel_position),
                                       erosion_py_width)
                # these aren't defined if pys evolve
        elif not evolving_pys:
            stepwise_depo = np.ones(nt, dtype=bool)
            stepwise_inc = np.ones(nt, dtype=bool)

        if restricted_channel_mass_conserved or decoupled_erosion_depo:
            depo_term_in = float(depo_py)  # attempt to force copy
            inc_term_in = float(erosion_py)
        elif walking_erosion_depo or compensation:
            depo_term_in = float(depo_py_width)
            inc_term_in = float(erosion_py_width)
        elif never_erosion:
            depo_term_in = 1.
            inc_term_in = 0.
        elif not modified_rules:
            depo_term_in = 1.
            inc_term_in = 1.

        R = 0.  # initial shoreline position
        doc = SL_trajectory[0]
        RF = doc/SF  # initial toe position
        Qvol = 0.  # sed vol added to delta at current time

        pr = 1 # variable for setting print interval

        diff_thresh = 0.001  # threshold for mass balancing
        relax_underest = 0.005

        tstore = []
        Rstore = []
        eroded_this_step = np.zeros(nt, dtype=bool)
        deposited_this_step = np.zeros(nt, dtype=bool)
        erosion_allowed_this_step = np.zeros(nt, dtype=bool)
        strat_eta = np.empty((nt, n-1), dtype=float)  # yikes!
        rollover_pixel = np.zeros(nt, dtype=int)
        rollover_eta = np.zeros_like(rollover_pixel, dtype=float)
        rollover_preserved = np.zeros((nt, nt), dtype=bool)
        self.rollover_preserved = rollover_preserved

        Qvol_lasttime = 0.
        DelVol = 0.

        for jj in xrange(nt):  # time loop
            print(jj)

            if compensation:
                pass

            if evolving_pys is True:
                # question arised here if we should also modify sedflux_modifier when
                # the channel becomes unconfined. I think... yes?
                #transiently update the pys for each case:
                if not modified_rules:
                    depo_can_occur = True
                    erosion_can_occur_on_topset = True
                elif never_erosion:
                    depo_can_occur = True
                    erosion_can_occur_on_topset = False
                elif restricted_channel_mass_conserved or decoupled_erosion_depo:
                    if using_rule_1:
                        if not channelised:
                            depo_py = 1.
                            erosion_py = 1.
                            sedflux_modifier = 1.
                            erosion_modifier = 1.
                        else:
                            depo_py = depo_term_in
                            erosion_py = inc_term_in
                            if decoupled_erosion_depo:
                                sedflux_modifier = 1./depo_term_in
                                erosion_modifier = inc_term_in/depo_term_in
                            else:
                                sedflux_modifier = 1./depo_term_in
                                erosion_modifier = 1.
                    if using_rule_2:
                        if not eroded_last_time:
                            depo_py = 1.
                            erosion_py = 1.
                            sedflux_modifier = 1.
                            erosion_modifier = 1.
                        else:
                            depo_py = depo_term_in
                            erosion_py = inc_term_in
                            if decoupled_erosion_depo:
                                sedflux_modifier = 1./depo_term_in
                                erosion_modifier = inc_term_in/depo_term_in
                            else:
                                sedflux_modifier = 1./depo_term_in
                                erosion_modifier = 1.
                    depo_can_occur = np.greater(py_thresholds_depo[jj],
                                                1.-depo_py)
                    if depo_can_occur:
                        erosion_can_occur_on_topset = np.greater(
                                py_thresholds_inc[jj], 1.-erosion_py/depo_py)
                    else:
                        erosion_can_occur_on_topset = False
                    # ...because we forbid erosion w/o depo
                elif walking_erosion_depo:
                    if using_rule_1:
                        if not channelised:
                            depo_py_width = 1.
                            erosion_py_width = 1.
                            sedflux_modifier = 1.
                            erosion_modifier = 1.
                        else:
                            depo_py_width = depo_term_in
                            erosion_py_width = inc_term_in
                            sedflux_modifier = 1./depo_term_in
                            erosion_modifier = inc_term_in/depo_term_in
                    if using_rule_2:
                        if not eroded_last_time:
                            depo_py_width = 1.
                            erosion_py_width = 1.
                            sedflux_modifier = 1.
                            erosion_modifier = 1.
                        else:
                            depo_py_width = depo_term_in
                            erosion_py_width = inc_term_in
                            sedflux_modifier = 1./depo_term_in
                            erosion_modifier = inc_term_in/depo_term_in
                    depo_can_occur = np.less(np.absolute(section_position -
                                                         channel_position[jj]),
                                             depo_py_width)
                    erosion_can_occur_on_topset = np.less(np.absolute(
                                                section_position-channel_position[jj]),
                                                          erosion_py_width)
                else:
                    raise NameError
            else:
                depo_can_occur = stepwise_depo[jj]
                erosion_can_occur_on_topset = stepwise_inc[jj]

            #Qvol += delt  # vol of sed added to delta  # DEJH added the Q
            if depo_can_occur or compensation:
                Qvol += Q*delt*sedflux_modifier
                R += delr.min()
                diff = 1.
                # ^Diff betw vol of sed added & vol in delta. The value of 1 is a
                # default to force the while loop to start. BUT only if depo allowed!
            else:
                Qerode = 0.
                diff = 0.


            doc = SL_trajectory[jj]  # set SL
            #R += delr  # initial guess for shoreline
            RF = R + doc/SF  # initial foreset toe position

            last_diff = -np.inf
            while abs(diff)> diff_thresh:
                # calculate new sed height in each volume using geometry rules

                distal_nodes = np.greater(rnode, RF)
                foreset_nodes = np.logical_and(np.less_equal(rnode, RF),
                                               np.greater(rnode, R))
                topset_nodes = np.less_equal(rnode, R)

                # in the foreset, sed height is set by geometry of foreset OR
                # the old height, whichever is largest
                # allows abandonment of the foreset
                # erosion can never occur

                etanew[distal_nodes] = eta[distal_nodes]  # ...always
                if depo_can_occur:
                    etanew[foreset_nodes] = np.maximum((RF-rnode[foreset_nodes])*SF,
                                                        eta[foreset_nodes])
                else:
                    etanew[foreset_nodes] = eta[foreset_nodes]

                # in the topset, sed height always by geometric rule if erosion
                # allowed, same max approach as above if it isn't
                possible_topset_eta = (R-rnode[topset_nodes])*ST + doc
                # ^ this is the height of the topset we would get if no other
                # strata present
                Qerode = np.sum((eta[topset_nodes]-possible_topset_eta
                    ).clip(0.)*delr[topset_nodes]*rnode[topset_nodes
                        ])*theta*inc_term_in
                Qin_tot = (Qerode + Q*delt)/(theta * depo_term_in)
                etanew[topset_nodes] = np.maximum(possible_topset_eta,
                                                  eta[topset_nodes])
                # ^ no erosion permitted yet...

                # accumulate the sed vol
                # essentially the product of col_area*r*theta
                # NOTE this will pick up previously abandoned foreset seds
                DelVol = (etanew*delr*rnode).sum()
                diff = Qvol_lasttime + Qin_tot - DelVol
                if erosion_can_occur_on_topset:
                    if depo_can_occur:
                        etanew[topset_nodes] = possible_topset_eta
                        DelVol = (etanew*delr*rnode).sum()
                    else:
                        raise NameError('This condition is forbidden!!')
                        # etanew[topset_nodes] = eta[topset_nodes]
                        # Qerode = 0.
                else:  # erosion forbidden
                    pass #this is already set right

                # calc the diff betw vol added & vol in delta
                # should be zero - a mismatch requires updating guess for shoreline


                # update guess for shoreline and toe
                # very crude - just driven by val of diff
                # use a constant for under-relaxation found by trial&error
                R += relax_underest*diff
                RF = R + doc/SF
                #print diff

            # loop terminates when diff is v small, indicating vol entered is same
            # as vol deposited

            Qvol_lasttime = float(DelVol)  # force a copy
            tstore.append(jj*delt)
            Rstore.append(R)

            # handle the stratigraphic preservation
            # This approach is memory hungry, but easy
            current_head = float(etanew[0])
            eta_at_R = np.interp(R, rnode, etanew)
            if jj != 0:
                prev_head = eta[0]
            else:
                prev_head = None
            if current_head >= prev_head:
                # this condition means a new topset starting at the wall
                # ...& downcutting did not occur
                downcut_this_time = False

            else:
                # either we downcut, or no erosion occurred and the topset is far out,
                # or no depo occurred either. Latter is a special case, where we
                # mustn't store anything at all this time.
                if erosion_can_occur_on_topset:
                    downcut_this_time = True
                else:
                    downcut_this_time = False

            if downcut_this_time:
                # Search and update the strat
                if print_strata:
                    strat_eta[:jj, :] = np.where(strat_eta[:jj, :] > etanew,
                                                 etanew, strat_eta[:jj, :])
                current_eta_old_rollovers = etanew[rollover_pixel[:jj]]
                #rollover_eta = np.where(current_eta_old_rollovers>rollover_eta[:jj],
                #                        rollover_eta, current_eta_old_rollovers)
                erased_rollovers = np.where(current_eta_old_rollovers<rollover_eta[:jj])
                #rollover_eta[:jj][erased_rollovers] = current_eta_old_rollovers[erased_rollovers]
                rollover_preserved[jj-1, erased_rollovers] = False
                #print 'ceor: ', current_eta_old_rollovers
                #print 'roll_eta: ', rollover_eta
                #print 'erased: ', erased_rollovers
            if print_strata:
                strat_eta[jj, :] = etanew.copy()
            rollover_preserved[jj, :] = rollover_preserved[jj-1, :].copy()
            if depo_can_occur:
                if R>=0.:  # R<0 => degenerate topo, no rollover
                    # could also argue there *is* preservation, but choose this way
                    #print R
                    rollover_preserved[jj, jj] = True
            rollover_pixel[jj] = np.searchsorted(px_left_edge, R)
            rollover_eta[jj] = etanew[rollover_pixel[jj]].copy()

            # end of tstep - switch old and new col sed heights
            # store position time and postion of shoreline
            eta = etanew.copy()
            # more bookkeeping
            if erosion_can_occur_on_topset:
                eroded_last_time = downcut_this_time
            erosion_allowed_this_step[jj] = erosion_can_occur_on_topset
            eroded_this_step[jj] = downcut_this_time
            deposited_this_step[jj] = depo_can_occur
            if evolving_pys and using_rule_1:
                if depo_can_occur:  # make no changes if no depo
                    axis_intercept = eta_at_R+ST*R
                    channelised = bool(np.argmax((axis_intercept, historic_high)))
                    historic_high = np.amax((axis_intercept, historic_high))

            if print_strata:
                close('all')
                #figure(3)
                #for i in xrange(nt):
                #    plot(rnode, strat_eta[i,:])
                #show()

        (timescales_of_completeness, completeness_at_multiple) = self.full_completeness()

        completeness_records.append(completeness_at_multiple)

        if graphs:
            #print out shoreline
            figure(2)
            plot(tstore, Rstore)

            figure(6)
            plot(tstore, SL_trajectory)

            #print the strat
            if print_strata:
                figure(3)
                for i in xrange(nt):
                    plot(rnode, strat_eta[i,:])

            figure(4)
            plot(rnode, eta)

            figure(5)
            plt.gca().set_xscale('log')
            plt.xlabel('multiple of smallest tstep')
            plt.ylabel('completeness')
            for i in completeness_records:
                plot(timescales_of_completeness, i)

        return timescales_of_completeness, completeness_records



    def read_input_file(self, input_file):
        d = {}
        assert type(input_file) is str
        with open(input_file) as f:
            for line in f:
               (key, val) = line.split()
               d[str(key)] = val
        return d

    def completeness_at_tscale(self, timescale_multiple=1):
        num_windows = self.nt-timescale_multiple+1
        completeness_count = 0
        for i in xrange(num_windows):
            completeness_count += np.any(self.rollover_preserved[-1,
                                                     i:(i+timescale_multiple)])
        return float(completeness_count)/num_windows

    def full_completeness(self, timestep=1.):
        """
        Calculates, stores and plots the full range of possible completenesses
        for the whole delta.
        """
        completeness_at_multiple = np.empty(self.nt-1, dtype=float)
        timescales_of_completeness = (np.arange(self.nt-1, dtype=float)+1
                                      )*self.delt
        for i in xrange(1, self.nt):
            completeness_at_multiple[i-1] = self.completeness_at_tscale(i)
        return timescales_of_completeness, completeness_at_multiple


if __name__ == "__main__":
    mydelta = delta()
    ins = mydelta.read_input_file('test_inputs.txt')
    nt = int(ins['nt'])
    SL_trajectory = (np.random.rand(nt)*2.+0.2)
    SL_trajectory += np.arange(nt)*0.01
    #SL_trajectory = np.ones(nt)
    mydelta.execute('test_inputs.txt', SL_trajectory)
    # completenesses = []
    # for i in xrange(10):
    #     SL_trajectory = (np.random.rand(nt)*2.+0.2)
    #     SL_trajectory += np.arange(nt)*0.01
    #     tscales, completenesses = mydelta.execute('test_inputs.txt', SL_trajectory, completeness_records=completenesses, graphs=False)
    # figure(1)
    # plt.gca().set_xscale('log')
    # plt.xlabel('multiple of smallest tstep')
    # plt.ylabel('completeness')
    # for i in completenesses:
    #     plot(tscales, i)
    show()
