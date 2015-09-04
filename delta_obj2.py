## This is a object-based version o Man Liang & Vaughan Voller's discritized
## delta xsectional model, after man_delta4.py.
## DEJ Hobley, Aug 2015.

# assumptions:
# The topset always maintains its slope - eroding if required
# In a rising BL the foreset can be abandoned
# erosion w/o depo is always forbidden in a single section

# Most likely cause of this code crashing is it reaching the RHS of its domain

import numpy as np
from pylab import figure, plot, show, close
import matplotlib.pyplot as plt

class delta(object):
    """
    This object creates, updates, and maintains a geometrically controlled
    synthetic delta topography.
    """

    init_inputs = set(['n', 'delr', 'delt', 'nt', 'Q', 'ST', 'SF', 'theta'])
    comp_inputs = set(['compensation_depth', 'erosion_py_width',
                       'depo_py_width', 'drift'])
    evolving_inputs = set(['rule',])
    restricted_inputs = set(['activity_py',])
    decoupled_inputs = set(['erosion_py', 'depo_py'])
    walking_inputs = set(['erosion_py_width', 'depo_py_width', 'drift'])

    def execute(self, input_file, SL_trajectory,
                print_strata=True, graphs=True, completeness_records=[],
                initial_topo=None,
                evolving_pys = False,
                restricted_channel_mass_conserved = False,
                decoupled_erosion_depo = False,
                walking_erosion_depo=False,
                compensation=False,
                never_erosion=False,
                **kwds):
        """
        Executes the construction of the geometric delta, following mass
        balance scheme modifications as specified by the arguments. The delta
        is represented as a one dimensional array, with the column height
        representing the topography in this virtual delta section. Node spacing
        is quadratic, to honor the spreading of the delta wedge in 2D.

        This delta is used to explore trends in stratigraphic completeness,
        and in particular, in shoreline completeness, sensu Mahon et al (2015).

        Baselevel is the tank base, assumed at height = 0.

        Parameters
        ----------
        input_file : string
                Specifies the path to the text input file containing the
                parameters with which to initialize the run.
                Compulsory parameters are
                    n : the number of nodes to use in the grid
                    delr : spacing of the first two nodes in the array
                    delt : the timestep at which to record stratigraphy
                    nt : the total number of timesteps
                    Q : the sediment flux, volume per time
                    ST : the slope of the topset
                    SF : the slope of the foreset
                    theta : the spreading angle of the delta wedge
                Additional parameters needed by each rule for building the
                delta are listed below.
        SL_trajectory : array-like
                the sea level at each timestep
        print_strata : bool (default True)
                Whether to record and plot the actual preserved
                stratigraphy in the delta. Setting False can accelerate
                the code.
        graphs : bool (default True)
                Whether to plot graphs visualizing output once the method
                completes. The plots are:
                1. sea level vs time
                2. shoreline position vs time
                3. delta cross sectional stratigraphy (unless suppressed)
                4. delta final topography
                5. shoreline completeness vs timescale diagram
                6. if 'compensation', a plot of virtual depositional
                   thickness around the arc of the delta (used to move the
                   channel line around to allow compensatory stacking)
        completeness_records : list (optional)
                If set, the new completeness array is appended to this list,
                and this list returned, rather than just the array being
                returned in a new, len-1 list.
        initial_topo : len-2 tuple, or len-n array (optional)
                If specified, allows an initial topography to be set. If an
                array, this array is the initial topography. If a len-2 tuple,
                this is the (r, z) coordinates of a rollover, defining a simple
                initial delta topography.
        evolving_pys : bool (default False)
                If True, implements one of two rules which make incisional
                dynamics different from aggradational dynamics. If rule==1,
                the historic maximum topography is tracked. If current topo
                is above this, deposition occurs across the whole delta top. If
                below, if is restricted to a narrow channel (controlled by the
                depo_py). If rule==2, the behavior in only the last timestep
                determines whether the whole delta top experiences deposition.

        These parameters are all then mutually exclusive. If none is set,
        deposition and erosion both occur uniformly over the whole delta top at
        all timesteps.

        restricted_channel_mass_conserved : bool (default False)
                If True, the parameter *activity_py* sets the probability that
                incision and erosion can occur at any given timestep. Each
                timestep is independent. If erosion can occur, deposition can
                also occur, and vice versa.
        decoupled_erosion_depo : bool (default False)
                If True, the erosional probability and depositional probability
                can be different. However, we still enforce that if erosion
                can occur, then deposition must also be allowed in the section
                (i.e., depo_py >= erosion_py). The probabilities are set by
                the parameters *depo_py* and *erosion_py*. Each timestep is
                independent.
        walking_erosion_depo : bool (default False)
                If True, the location of the channel center line, about which
                the erosion and deposition loci are assumed to be symmetrically
                distributed, is allowed to perform a random walk over the delta
                top. The erosion and deposition widths are set with
                *depo_py_width* and *erosion_py_width*, and both are scaled to
                1 (i.e., they remain probabilities). The random walk takes
                uniformly randomly distributed steps up to a value of *drift*
                (also scaled to 1). Boundary conditions are looped.
        compensation : bool (default False)
                If true, normal behavior is like walking_erosion_depo, and the
                same parameters must be set. However, it also takes the param
                *compensation_depth*, which allows "avulsions" on the delta
                top. This is executed by tracking building a ghost topography
                beneath the channel center line; once this topography has
                aggraded more than compensation_depth since the last avulsion
                event, the channel center line moves. The switch occurs to a
                place on the delta top which is lowest or joint lowest (this
                is recorded by tracking a radial section around the delta, and
                assuming unit deposition of sediment thickness occurs
                symmetrically distributed around the channel line each time the
                channel avulses). After avulsion, the ghost section is assumed
                to now have initial topography equivalent to the real section
                at that time. This is broadly equivalent to the methods used
                by Straub et al. 2009 (JSR).
        never_erosion : bool (default False)
                If True, erosion can never occur, and deposition can always
                occur.

        Finally, note any parameter that must be set in the parameter file can
        also be passed as an argument to execute. This is intended to make
        parameter space explorations much easier.

        Returns
        -------
        timescales_of_completeness : array
                The multiples of delt represented by each entry of the arrays
                stored in completeness_records.
        completeness_records : list of arrays
                If the argument completeness_records was not supplied to
                execute(), this is a length-1 list containing an array giving
                the shoreline completeness at each timescale specified in
                timescales_of_completeness. If completeness_records was set,
                this input list will be returned, appended with this new array.


        """
        # load the input params from the textfile:
        input_dict = self.read_input_file(input_file)
        # also check we don't need to overwrite anything from the arguments...
        for varname in self.init_inputs:
            try:
                exec(varname + " = float(input_dict['" + varname + "'])")
            except KeyError:
                exec(varname + " = float(kwds['" + varname + "'])")
            try:  # if it's specified as a kwd, that takes precedence
                exec(varname + " = float(kwds['" + varname + "'])")
            except KeyError:
                pass
        # set the ones that should be ints:
        n = int(n)
        nt = int(nt)

        self.delt = delt
        self.nt = nt

        rnode = np.arange(0.5, n+0.5, 1.)**0.5 * delr/np.sqrt(0.5)
        px_left_edge = np.zeros_like(rnode)
        px_left_edge[1:] = (rnode[1:] - rnode[:-1])/2. + rnode[1:]
        delr = rnode[1:]-rnode[:-1] #this is cludgey, should do accurately
        rnode = rnode[:-1]
        sedflux_modifier = 1.
        erosion_modifier = 1.

        assert np.any(SL_trajectory<=0.) == False

        # check for initial topo:
        if initial_topo is None:
            etanew = np.zeros_like(rnode)  # sed height in column
            eta = np.zeros_like(rnode)
        elif len(initial_topo) == 2:  # (x,z) of a single rollover
            doc = initial_topo[1]
            R = initial_topo[0]
            RF = R + doc/SF
            distal_nodes = np.greater(rnode, RF)
            foreset_nodes = np.logical_and(np.less_equal(rnode, RF),
                                           np.greater(rnode, R))
            topset_nodes = np.less_equal(rnode, R)
            eta = np.zeros_like(rnode)
            eta[foreset_nodes] = (RF-rnode[foreset_nodes])*SF
            eta[topset_nodes] = (R-rnode[topset_nodes])*ST + doc
            etanew = eta.copy()
        else:  # a whole topo was supplied
            assert len(initial_topo) == rnode.size
            eta = np.array(initial_topo)
            etanew = eta.copy()

        # now the stuff conditional on the param rules:
        if compensation:
            for varname in self.comp_inputs:
                try:
                    exec(varname + " = float(input_dict['" + varname + "'])")
                except KeyError:
                    exec(varname + " = float(kwds['" + varname + "'])")
                try:  # if it's specified as a kwd, that takes precedence
                    exec(varname + " = float(kwds['" + varname + "'])")
                except KeyError:
                    pass
            assert erosion_py_width <= depo_py_width
            sedflux_modifier /= depo_py_width
            erosion_modifier = erosion_py_width/depo_py_width  # as above
            # ^these need to get halved for comparison with pys above
            pos_or_neg = 1.  #make this -1. to reverse the drift direction
            channel_position_steps = 2.*(np.random.rand(nt)-0.5)*drift  # the steps
            channel_start = np.random.rand()
            channel_position = channel_position_steps.copy()
            channel_position[0] = channel_start
            section_position = np.random.rand()

            # ^ this tells us if we're counting the init pos as "active" or
            # "inactive" depo, and will govern how and when we jump the depo
            # locus.
            # \/ this stuff needs to go into the active loop now
            looped_round = np.zeros(nt, dtype=int)
            for i in xrange(nt-1):
                # a looped model makes most sense here
                ##channel_position[i+1] += pos_or_neg*channel_position[i]
                channel_position[i+1] += channel_position[i]
                #excess = channel_position[i+1]%1.
                if channel_position[i+1] > 1.:  # loop to the other end
                    ##pos_or_neg *= -1.
                    ##channel_position[i+1] = 1. - excess
                    channel_position[i+1] -= 1.
                    looped_round[i+1] = 1
                elif channel_position[i+1] < 0.:  # reverse direction from 0.
                    ##pos_or_neg *= -1.
                    ##channel_position[i+1] = 1. - excess
                    channel_position[i+1] += 1.
                    looped_round[i+1] = -1.
            compensational_base = 0.
            accumulated_OOS_py = 0.  # this is a ratchet to control when to
            # skip in & out of section
        if evolving_pys:
            assert compensation is False, "Cannot use evolving_pys with compensation!!"
            if not rule_override:
                for varname in self.evolving_inputs:
                    try:
                        exec(varname + " = float(input_dict['" + varname + "'])")
                    except KeyError:
                        exec(varname + " = float(kwds['" + varname + "'])")
                    try:  # if it's specified as a kwd, that takes precedence
                        exec(varname + " = float(kwds['" + varname + "'])")
                    except KeyError:
                        pass
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
            for varname in self.restricted_inputs:
                try:
                    exec(varname + " = float(input_dict['" + varname + "'])")
                except KeyError:
                    exec(varname + " = float(kwds['" + varname + "'])")
                try:  # if it's specified as a kwd, that takes precedence
                    exec(varname + " = float(kwds['" + varname + "'])")
                except KeyError:
                    pass
            erosion_py = activity_py
            depo_py = activity_py
            sedflux_modifier /= depo_py
        if decoupled_erosion_depo:
            for varname in self.decoupled_inputs:
                try:
                    exec(varname + " = float(input_dict['" + varname + "'])")
                except KeyError:
                    exec(varname + " = float(kwds['" + varname + "'])")
                try:  # if it's specified as a kwd, that takes precedence
                    exec(varname + " = float(kwds['" + varname + "'])")
                except KeyError:
                    pass
            assert erosion_py <= depo_py
            # ^we need to assume this to handle mass balance
            sedflux_modifier /= depo_py
            erosion_modifier = erosion_py/depo_py
        if walking_erosion_depo:
            for varname in self.walking_inputs:
                try:
                    exec(varname + " = float(input_dict['" + varname + "'])")
                except KeyError:
                    exec(varname + " = float(kwds['" + varname + "'])")
                try:  # if it's specified as a kwd, that takes precedence
                    exec(varname + " = float(kwds['" + varname + "'])")
                except KeyError:
                    pass
            assert erosion_py_width <= depo_py_width
            sedflux_modifier /= depo_py_width
            erosion_modifier = erosion_py_width/depo_py_width  # as above
            section_position = np.random.rand()
            # ^these need to get halved for comparison with pys above
            assert drift <= 1.
            ##pos_or_neg = 1.  #make this -1. to reverse the drift direction
            channel_position_steps = 2.*(np.random.rand(nt)-0.5)*drift  # the steps
            channel_start = np.random.rand()
            channel_position = channel_position_steps.copy()
            channel_position[0] = channel_start
            looped_round = np.zeros(nt, dtype=int)
            for i in xrange(nt-1):
                # a looped model makes most sense here
                ##channel_position[i+1] += pos_or_neg*channel_position[i]
                channel_position[i+1] += channel_position[i]
                #excess = channel_position[i+1]%1.
                if channel_position[i+1] > 1.:  # loop to the other end
                    ##pos_or_neg *= -1.
                    ##channel_position[i+1] = 1. - excess
                    channel_position[i+1] -= 1.
                    looped_round[i+1] = 1
                elif channel_position[i+1] < 0.:  # reverse direction from 0.
                    ##pos_or_neg *= -1.
                    ##channel_position[i+1] = 1. - excess
                    channel_position[i+1] += 1.
                    looped_round[i+1] = -1.
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
                channel_pos_off_ends = np.where(looped_round!=0,
                                                channel_position+looped_round,
                                                channel_position)
                stepwise_depo = np.less(np.absolute(section_position-channel_pos_off_ends),
                                        depo_py_width/2.)
                stepwise_inc = np.less(np.absolute(section_position-channel_pos_off_ends),
                                       erosion_py_width/2.)
                if compensation:
                    deposition_accum = np.zeros(100, dtype=int)
                    # ^ a quantized number line to record where and how much
                    # depo we've had...
                    depo_in_section = stepwise_depo[0]
                    if depo_in_section:
                        accumulated_OOS_py_incl_section = True
                    else:
                        accumulated_OOS_py_incl_section = False
                    channel_position_at_switch = channel_start
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
        R_OOS = 0.  # init for the alt. section, if using compensation
        doc = SL_trajectory[0]
        RF = doc/SF  # initial toe position
        RF_OOS = float(RF)
        if initial_topo is None:
            Qvol = 0.  # sed vol added to delta at current time
        else:
            Qvol = (eta*delr*rnode).sum()

        pr = 1 # variable for setting print interval

        diff_thresh = 0.0005*Q**0.333/delt  # threshold for mass balancing
        relax_underest = 0.01/delr[0]/Q**0.333/delt

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

        Qvol_lasttime = float(Qvol)
        DelVol = 0.

        if compensation:
            Qvol_lasttime_OOS = float(Qvol)
            out_of_section_eta = eta.copy()
            out_of_section_oldeta = eta.copy()
            OOS_head_baselevel = out_of_section_oldeta[0]  # updates below
            # the scheme here is that if depo_in_section is FALSE at t=0, we
            # track BOTH the actual section and also the topo at that point.
            # once *that* section hits the compensation depth, we "avulse" to
            # a new randomly selected section, resetting the topo to whatever
            # it was at the time of the first avulsion. Once we have "filled"
            # the "other topo", we avulse to a low point.
            # This is achieved by tracking in some sense where depo is
            # happening between avulsion events, using the deposition_accum
            # array, then selecting the lowest areas from this.
            # Note we move the section of we're tracking *with* the stream
            # line. i.e., the compensation depth doesn't feel the "walking"
            # of the channel position, and it always goes up and down if it
            # can.

        for jj in xrange(nt):  # time loop
            print(jj)

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
                                                     channel_pos_off_ends[jj]),
                                             depo_py_width/2.)
                    erosion_can_occur_on_topset = np.less(np.absolute(
                                    section_position-channel_pos_off_ends[jj]),
                                                          erosion_py_width/2.)
                else:
                    raise NameError
            else:
                depo_can_occur = stepwise_depo[jj]
                erosion_can_occur_on_topset = stepwise_inc[jj]

            #Qvol += delt  # vol of sed added to delta  # DEJH added the Q
            if depo_can_occur:
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

            # now a special case for if we're doing compensation.
            # the above still applies - in the actual section - but we're also
            # handling a second section, which tracks the evolution of an
            # idealized "out of section" delta slice
            # Note that this section is treated as always eroding and always
            # depositing... Visualized as effectively moving the section line
            # with the channel thalweg as it migrates.
            if compensation:
                # Qvol is inherited from above
                R_OOS += delr.min()
                diff = 1.
                RF_OOS = R_OOS + doc/SF
                while abs(diff)> diff_thresh:
                    distal_nodes = np.greater(rnode, RF_OOS)
                    foreset_nodes = np.logical_and(np.less_equal(rnode, RF_OOS),
                                                   np.greater(rnode, R_OOS))
                    topset_nodes = np.less_equal(rnode, R_OOS)
                    out_of_section_eta[distal_nodes] = out_of_section_oldeta[distal_nodes]
                    out_of_section_eta[foreset_nodes] = np.maximum((RF_OOS-rnode[foreset_nodes])*SF,
                                                        out_of_section_oldeta[foreset_nodes])
                    # ^...depo now always occurs
                    possible_topset_OOS_eta = (R_OOS-rnode[topset_nodes])*ST + doc
                    Qerode_OOS = np.sum((out_of_section_oldeta[topset_nodes]-possible_topset_OOS_eta
                        ).clip(0.)*delr[topset_nodes]*rnode[topset_nodes
                            ])*theta*inc_term_in
                    Qin_tot_OOS = (Qerode_OOS + Q*delt)/(theta * depo_term_in)
                    out_of_section_eta[topset_nodes] = np.maximum(possible_topset_OOS_eta,
                                                      out_of_section_eta[topset_nodes])
                    DelVolOOS = (out_of_section_eta*delr*rnode).sum()
                    diff = Qvol_lasttime_OOS + Qin_tot_OOS - DelVolOOS
                    out_of_section_eta[topset_nodes] = possible_topset_OOS_eta
                    # DelVol = (out_of_section_eta*delr*rnode).sum()  # don't think now needed
                    R_OOS += relax_underest*diff
                    RF_OOS = R_OOS + doc/SF
                Qvol_lasttime_OOS = float(DelVolOOS)
                current_head_OOS = float(out_of_section_eta[0])
                # Now here, we modify the section position for next &
                # subsequent times if we exceed the compensation depth
                # (but the below all totally still works for the actual
                # section!):
                head_elevation = current_head_OOS - OOS_head_baselevel
                superelevated = head_elevation > compensation_depth
                if superelevated:
                    print "superelevated!"
                    # first, elevate the section we just worked on:
                    high_edge = int((channel_position_at_switch+depo_py_width/2.)*100)
                    low_edge = int((channel_position_at_switch-depo_py_width/2.)*100)
                    deposition_accum[np.clip(low_edge, 0, 100):np.clip(high_edge, 0, 100)] += 1
                    if high_edge >= 100:
                        assert high_edge < 200
                        deposition_accum[:(high_edge-100)] += 1
                    elif low_edge < 0:
                        assert low_edge >= -100
                        deposition_accum[(low_edge+100):] += 1
                    # then it's switchin' time
                    # find the possible minima places...
                    lowest_elev_val = deposition_accum.min()
                    poss_places = np.where(deposition_accum==lowest_elev_val)[0]
                    # pick a low one...
                    which_one = int(np.random.rand()*poss_places.size)
                    channel_position_at_switch = (poss_places[which_one]+np.random.rand())/100.
                    if jj != nt-1:
                        channel_position[jj+1] = float(channel_position_at_switch)
                        looped_round[jj+1] = 0
                        # reset the form to whatever the actual section
                        # currently looks like!
                        out_of_section_oldeta[:] = eta  # copy is automatic
                        OOS_head_baselevel = out_of_section_oldeta[0]
                else:
                    if jj!=nt-1:  # all as above to reset new ch positions
                        channel_position[jj+1] = channel_position[jj] + \
                                                   channel_position_steps[jj+1]
                        channel_pos_off_ends[jj+1] = float(channel_position[jj+1])
                        if channel_position[jj+1] > 1.:  # loop to the other end
                            channel_position[jj+1] -= 1.
                            looped_round[jj+1] = 1
                        elif channel_position[jj+1] < 0.:  # reverse direction from 0.
                            channel_position[jj+1] += 1.
                            looped_round[jj+1] = -1.
                        stepwise_depo[jj+1] = np.less(np.absolute(section_position -
                                                  channel_pos_off_ends[jj+1]),
                                                  depo_py_width/2.)
                        stepwise_inc[jj+1] = np.less(np.absolute(section_position -
                                                  channel_pos_off_ends[jj+1]),
                                                  erosion_py_width/2.)
                        out_of_section_oldeta[:] = out_of_section_eta

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
                erased_rollovers = np.where(current_eta_old_rollovers<rollover_eta[:jj])[0]
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

        self._volume = DelVol*theta
        self.tstore = tstore
        self.Rstore = Rstore
        self.strat_eta = strat_eta
        self.rnode = rnode

        if graphs:
            #SL
            figure(1)
            plot(self.output_times, SL_trajectory)
            plot(self.output_times[self.final_preserved],
                 SL_trajectory[self.final_preserved], 'r.')
            plt.ylabel('Sea level')
            plt.xlabel('Time')

            #shoreline position
            figure(2)
            plot(self.output_times, self.shoreline_positions)
            plot(self.output_times[self.final_preserved],
                 self.shoreline_positions[self.final_preserved], 'r.')
            plt.ylabel('Shoreline position')
            plt.xlabel('Time')

            #print the strat
            if print_strata:
                figure(3)
                for i in xrange(nt):
                    plot(rnode, strat_eta[i,:])
                plt.xlabel('Radial distance')
                plt.ylabel('Height')

            #print the final topography
            figure(4)
            plot(rnode, eta)
            plt.xlabel('Radial distance')
            plt.ylabel('Topographic height')

            #print the completeness vs timecale diagram
            figure(5)
            plt.gca().set_xscale('log')
            plt.xlabel('multiple of smallest tstep')
            plt.ylabel('completeness')
            for i in completeness_records:
                plot(timescales_of_completeness, i)

            # if appropriate, plot the virtualized loci of deposition around
            # the delta arc (used in compensational stacking)
            if compensation:
                figure(6)
                plot(deposition_accum)
                plt.ylabel('depo thickness map')
                plt.xlabel('percent distance around delta arc')

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

    @property
    def volume(self):
        return self._volume

    @property
    def final_preserved(self):
        return self.rollover_preserved[-1]

    @property
    def output_times(self):
        return np.array(self.tstore)

    @property
    def shoreline_positions(self):
        return np.array(self.Rstore)

    @property
    def radial_distances(self):
        return self.rnode

    @property
    def strata(self):
        return self.strat_eta


if __name__ == "__main__":
    mydelta = delta()
    ins = mydelta.read_input_file('test_inputs.txt')
    nt = int(ins['nt'])
    SL_trajectory = (np.random.rand(nt)*2.+0.2)
    SL_trajectory += np.arange(nt)*0.01
    mydelta.execute('test_inputs.txt', SL_trajectory)
    show()
