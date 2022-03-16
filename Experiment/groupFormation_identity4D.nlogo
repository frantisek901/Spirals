;; Model for study of mutual interaction between public opinion and network structure.
;; RQ: How polarisation of public opinion shapes the network structure
;;     and how the network structure shapes public opinion and influences individual opinions?
;;
;; This code is from project of Mike, Ashley, Ashwin and FranCesko,
;; we apply Hegselmann-Krausse model in more than 1D and look how agents adapt in >1D opinion space and whether they form groups,
;; then we include small-world network (Watts-Strogatz) as another constraint, features of spiral of silence, and
;; individual assignment of uncertainity, tollerance, comformity and outspokeness.
;;
;; !!!EXPERIMENTAL FEATURE: IDENTITY!!!
;; !!!SUPER EXPERIMENTAL: INDIVIDUAL IDENTITY VISION!!!
;;

;; Created:  2021-10-21 FranCesko
;; Edited:   2022-03-16 FranCesko
;; Encoding: windows-1250
;; NetLogo:  6.2.2
;;

;; IDEA: What about simply employ Spiral of Silence?
;;       Just simply -- general parameter on scale (0; 1> and probability of speaking her attitude/opinion,
;;       baseline is p==1, everybody speaks always, if p==0.5 so everybody has 0.5 probability to speak her opinion/attitude at given step,
;;       if succeeds - speaks in given step, if not - falls silent for the respective step.
;;       In HK mechanism, agent computes mean opinion of all speaking agents who are inside 'opinion boundary' (are not further than threshold).
;;       In Defuant, agent randomly takes one speaking agent inside the 'opinion boundary' and sets opinon as average of their opinions.
;; DONE!
;;
;; IDEA: Handle P-speaking as Uncertainty -- besides constant value for every agent, create random mode (random uniform for the start),
;;       where all agents will have their own value of speaking probability which they will follow.
;; DONE!
;;
;; IDEA: Choose, how many opinions agents update: sometimes 1, 2, 3, 4 ...
;; DONE!
;;
;; IDEA: Employ Schelling principle -- if the agents are unhappy in their neighborhood,
;;       the cut off all the links and create new set of links, i.e., join new neighborhood.
;; DONE!
;;
;; IDEA: Give weights to opinions... Taken from media, or from interpersonal communication:
;;          - agents pick opinion according the importance, and update importance according number of contacts regarding the opinion
;;
;; IDEA: Compute clusters.
;; Elle: I do cluster detection using igraph:: cluster_walktrap() in R
;;


;; WISHLIST:
;;     - differentiate between interpersonal communication and social media communication -- two overlapping networks with their own rules
;;     - how radicalization is possible? How polarization happens?
;;     - differential attraction
;;     - repulsion
;;     - media exposure will be crucial…we can ask abt opinion consistent content, opinion contrary, and “mainstream/mixed”…
;;       how to we conceptualize/model those in ABM? Is this too simplistic (eg, think of the different flavors of conservative media,
;;       ranging from CDU type media to extremist hate groups).
;;     - how to think about social media influencers (eg Trump before deplatforming)…
;;       is it possible to designate “superagents” who influence everyone sharing certain beliefs and see their effects…
;;       both reach everyone in a group and their opinions are very highly weighted (or people vary in how much they weight that opinion?
;;       Could estimate Twitter effect that way!  Perhaps one could even model how movement towards an opinion might influence the superagent
;;       to increase communication or change focus…
;;     - Employ homophily/heterophily principle at model start.
;;     - Control degree of opinion randomness at the start (different mean and SD of opinion for different groups)
;;     - Mike was thinking…after we do “superagents”, the Trump/foxnews avatars…one thing that would be neat and represent social reality
;;       is to have some kind of attraction to those who share beliefs (including superagents), but that decreases with close proximity…
;;       that way we have less ability/willingness to select attitude consistent sources around us (eg can’t escape family and coworkers),
;;       but can seek them elsewhere.  That might allow us to look at what happens in a more or less diverse local opinion environment, among other things.
;;     - Use clustering algorithm for creating group identities


;; TO-DO:
;; 1) constructing file name for recording initial and final state of simulation
;; DONE!
;; 2) implementing recording into the model -- into setup and final steps (delete component detection and just record instead)
;; DONE!
;;
;; 3) Reviewer's comments:
;; The reviewer for your computational model Simulating Components of the Reinforcing Spirals Model and Spiral of Silence v1.0.0 has recommended that changes be made to your model. These are summarized below:
;; Very interesting model! It needs better documentation though, both within the code as comments, and the accompanying narrative documentation. Please consider following the ODD protocol or equivalent to describe your model in sufficient detail so that another could replicate the model based on the documentation.
;; Has Clean Code:
;; The code should be cleaned up and have more comments describing the intent and semantics of the variables.
;; Has Narrative Documentation:
;; The info tab is empty and the supplementary doc does not include sufficient detail to replicate the model. For example documentation please see ODD examples from other peer reviewed models in the library.
;; Is Runnable:
;; Runs well.
;; On behalf of the editors@comses.net, thank you for submitting your computational model(s) to CoMSES Net! Our peer review service is intended to serve the community and we hope that you find the requested changes will improve your model’s accessibility and potential for reuse. If you have any questions or concerns about this process, please feel free to contact us.
;;
;; 4) Adapt recording data for cluster computation -- machine's root independent.
;; DONE!
;;
;; 5) Appropriate recorded data format -- we want it now as:
;;    a) dynamical multilayer network, one row is one edge of opinion distance network,
;;    b) separate file with agent's traits (P-speaking, Uncertainty etc.)
;;    c) as it was before, contextual variables of one whole simulation run are coded in the filenames
;; DONE!
;;
;; 6) Implement K-clusters algorithm for addressing just 2 clusters.
;;
;; 7) Store state of simulation in global variables -- instead of storing it to file,
;;    we will to store states in variables and then at the end let the experiment/BehaviorSpace itself to record it to the file

extensions [nw]

breed [centroids centroid]
undirected-link-breed [comms comm]
undirected-link-breed [l-distances l-distance]

turtles-own [Opinion-position P-speaking Speak? Uncertainty Record Last-opinion
  Tolerance Conformity Satisfied? group distance_to_centroid Group-threshold Identity-group Opponents-ratio id-threshold-level]
l-distances-own [l-weight]
comms-own [op-weight]

globals [main-Record components positions network-changes agents positions_clusters
  polarisation normalized_polarisation unweighted_polarisation unweighted_normalized_polarisation ESBG_polarisation id_threshold_set
  betweenness_start eigenvector_start clustering_start modularity_start mean_path_start normalized_polarization_start ESBSG_polarization_start
  mean_op1_start mean_op2_start sd_op1_start sd_op2_start median_op1_start median_op2_start
  lower_op1_start lower_op2_start upper_op1_start upper_op2_start
  mean_op3_start mean_op4_start sd_op3_start sd_op4_start median_op3_start median_op4_start
  lower_op3_start lower_op4_start upper_op3_start upper_op4_start
  betweenness_final eigenvector_final clustering_final modularity_final mean_path_final normalized_polarization_final ESBSG_polarization_final
  mean_op1_final mean_op2_final sd_op1_final sd_op2_final median_op1_final median_op2_final
  lower_op1_final lower_op2_final upper_op1_final upper_op2_final
  mean_op3_final mean_op4_final sd_op3_final sd_op4_final median_op3_final median_op4_final
  lower_op3_final lower_op4_final upper_op3_final upper_op4_final
]


;; Initialization and setup
to setup
  ;; Redundant conditions which should be avoided -- if the boundary is drawn as constant, then is completely same whether agents vaguely speak or openly listen,
  ;; seme case is for probability speaking of 100%, then it's same whether individual probability is drawn as constant or uniform, result is still same: all agents has probability 100%.
  if avoid-redundancies? and mode = "vaguely-speak" and boundary-drawn = "constant" [stop]
  if avoid-redundancies? and p-speaking-level = 1 and p-speaking-drawn = "uniform" [stop]
  ;; these two conditions cover 7/16 of all simulations, approx. the half! This code should stop them from running.

  ;; We erase the world and clean patches
  ca
  ask patches [set pcolor patch-color]

  ;; We initialize small-world network with random seed
  if set-seed? [random-seed RS]
  if HK-benchmark? [set n-neis (N-agents - 1) / 2]
  nw:generate-watts-strogatz turtles comms N-agents n-neis p-random [
    fd (max-pxcor - 1)
    set size (max-pxcor / 10)
  ]

  ;; To avoid some random artificialities due to small-world network generation,
  ;; we have to set random seed again.
  if set-seed? [random-seed RS]
  ;; Then we migh initialize agents/turtles
  ask turtles [
    set Opinion-position n-values opinions [precision (1 - random-float 2) 3]  ;; We set opinions...
    set Last-opinion Opinion-position  ;; ...set last opinion as present opinion...
    set Record n-values record-length [0]  ;; ... we prepare indicator of turtle's stability, at all positions we set 0 as non-stability...
    set P-speaking get-speaking  ;; ...assigning individual probability of speaking...
    set speak? speaking  ;; ...checking whether agent speaks...
    set Uncertainty get-uncertainty  ;;... setting value of Uncertainty.
    set Tolerance get-tolerance  ;; Setting individual tolerance level, as well as ...
    set Conformity get-conformity  ;; setting individual conformity level, and ...
    set Group-threshold get-group-threshold  ;; Individual sensitivity for group tightness/threshold.
    set Identity-group no-turtles  ;; Now, we just initialize it, later may be... we use it also here meaningfully.
    set Opponents-ratio 0  ;; Fraction of opponents
    getColor  ;; Coloring the agents according their opinion.
    getPlace  ;; Moving agents to the opinion space according their opinions.
  ]
  set agents turtle-set turtles  ;; Note: If we just write 'set agents turtles', then variable 'agents' is a synonym for 'turtles', so it will contain in the future created centroids!
  ask agents [create-l-distances-with other agents]  ;; Creating full network for computing groups and polarisation
  ask l-distances [set hidden? true]  ;; Hiding links for saving comp. resources

  compute-identity-thresholds
  ;; Coloring patches according the number of agents/turtles on them.
  ask patches [set pcolor patch-color]
  ;; Hiding links so to improve simulation speed performance.
  ask comms [set hidden? TRUE]
  ;; Setting the indicator of change for the whole simulation, again as non-stable.
  set main-Record n-values record-length [0]
  ;; Setting communication and distances links' weights
  update-links-weights
  ;; Setting agents' identity groups
  ifelse use_identity? [
    if identity_type = "global" [set-group-identities]
    if identity_type = "individual" [set-individual-group-identities]
  ][
    ask agents [set Identity-group agents]
  ]
  ;; update satisfaction
  ask agents [set Satisfied? get-satisfaction]

  ;; Setting control variable of network changes
  set network-changes 0
  ;; Compute polarisation
  ;compute-polarisation-repeatedly  ;; We are computing it in next command after reseting ticks, so we save computer time here.

  reset-ticks

  ;;;; Finally, we record initial state of simulation
  ;; If we want we could construct filename to contain all important parameters shaping initial condition, so the name is unique stamp of initial state!
  if construct-name? [set file-name-core (word RS "_" N-agents "_" p-random "_" n-neis "_" opinions "_" updating "_" boundary "_" boundary-drawn "_" p-speaking-level "_"  p-speaking-drawn "_" mode)]
  ;; recording itself
  compute-initial-macro-state-of-simulation
end


;; Just envelope for updating agent at the begining of GO procedure
to set-group-identities
  ;; Cleaning environment
  ask centroids [die]

  ;; Detection of clusters via Louvain: Detection itself
  let selected-agents agents with [2 < count my-l-distances with [l-weight >= id_threshold]]  ;; Note: We take into account only not loosely connected agents
  nw:set-context selected-agents l-distances with [l-weight >= id_threshold]  ;; For starting centroids we take into account only not loosely connected agents, but later we set groups for all.
  let communities nw:louvain-communities
  ;repeat N-agents [set communities nw:louvain-communities]
  set N_centroids length communities

  ;; Computing clusters' mean 'opinion-position'
  set positions_clusters [] ;; List with all positions of all clusters
  foreach communities [c ->
    let one []  ;; List for one position of one cluster
    foreach range opinions [o -> set one lput precision (mean [item o opinion-position] of c) 3 one]
    set positions_clusters lput one positions_clusters
  ]

  ;; Preparation of centroids -- feedeing them with communities
  create-centroids N_centroids [
    set heading (who - min [who] of centroids)
    set Opinion-position item heading positions_clusters  ;; We set opinions, we try to do it smoothly...
    set shape "circle"
    set size 1.5
    set color 5 + (who - min [who] of centroids) * 10
    getPlace
  ]

  ;; Assignment of agents to groups
  ask agents [set group [who] of min-one-of centroids [opinion-distance]]  ;; Sic! Here we intentionally use all agents, including loosely connected.

  ;; Computation of centroids possitions
  compute-centroids-positions (agents)

  ;; Iterating cycle -- looking for good match of centroids
  while [sum [opinion-distance3 (Last-opinion) (Opinion-position)] of centroids > Centroids_change] [

    ;; turtles compute whether they are in right cluster and
    ask agents [set group [who] of min-one-of centroids [opinion-distance]]

    ;; Computation of centroids possitions
    compute-centroids-positions (agents)
  ]

  ;; Saving Identity group as agent-set
  ask agents [
    set Identity-group agents with [group = [group] of myself]
    if count Identity-group < 3 [set Identity-group agents]
  ]

  ;; Killing centroids without connected agents
  ask centroids [
    let wom who
    if (not any? agents with [group = wom]) [die]
  ]
  set N_centroids count centroids

  ;; Final coloring and killing of centroids
  if centroid_color? [ask agents [set color (5 + 10 * (group - min [who] of centroids))]]
  if killing_centroids? [ask centroids [die]]
end




;to set-individual-group-identities
;  ask agents [
;    ;; Cleaning environment
;    set Identity-group no-turtles
;    let my-Group-threshold Group-threshold
;
;    ;; Detection of clusters via Louvain: Detection itself
;    ;let selected-agents agents with [2 < count my-l-distances with [l-weight >= my-Group-threshold]]  ;; Note: We take into account only not loosely connected agents
;    nw:set-context agents l-distances with [l-weight >= my-Group-threshold]  ;; For starting centroids we take into account only not loosely connected agents, but later we set groups for all.
;    let communities nw:louvain-communities  ;; Louvain detection of communitites
;    foreach communities [c -> if member? self c [set Identity-group c]]  ;; Looking for 'self' in communities -- community which includes '(my)self' is set as Identity group.
;    if count Identity-group < 3 [set Identity-group agents]  ;; CHECK: If Identity group is (almost) empty, then set all agents as Identity group
;  ]
;end


; Setting identity groups individually via threshold levels to account for differences in sensitivity to group relationships
to set-individual-group-identities
  foreach range identity_levels[ i ->
    let idthresh item i id_threshold_set
    ;; Cleaning environment
    ask centroids [die]

    ;; Detection of clusters via Louvain: Detection itself
    let selected-agents agents with [2 < count my-l-distances with [l-weight >= idthresh]]  ;; Note: We take into account only not loosely connected agents
    nw:set-context selected-agents l-distances with [l-weight >= idthresh]  ;; For starting centroids we take into account only not loosely connected agents, but later we set groups for all.
    let communities nw:louvain-communities

    ;repeat N-agents [set communities nw:louvain-communities]
    set N_centroids length communities

    ;; Computing clusters' mean 'opinion-position'
    set positions_clusters [] ;; List with all positions of all clusters
    foreach communities [c ->
      let one []  ;; List for one position of one cluster
      foreach range opinions [o -> set one lput precision (mean [item o opinion-position] of c) 3 one]

      set positions_clusters lput one positions_clusters
    ]


    ;; Preparation of centroids -- feedeing them with communities
    create-centroids N_centroids [
      set heading (who - min [who] of centroids)
      set Opinion-position item heading positions_clusters  ;; We set opinions, we try to do it smoothly...
      set shape "circle"
      set size 1.5
      set color 5 + (who - min [who] of centroids) * 10
      getPlace
    ]

    ;; Assignment of agents to groups
    ask agents [set group [who] of min-one-of centroids [opinion-distance]]  ;; Sic! Here we intentionally use all agents, including loosely connected.


    ;; Computation of centroids possitions
    compute-centroids-positions (agents)

    ;; Iterating cycle -- looking for good match of centroids
    while [sum [opinion-distance3 (Last-opinion) (Opinion-position)] of centroids > Centroids_change] [

      ;; turtles compute whether they are in right cluster and
      ask agents [set group [who] of min-one-of centroids [opinion-distance]]

      ;; Computation of centroids possitions
      compute-centroids-positions (agents)
    ]

    ;; Saving Identity group as agent-set
    ask agents with [id-threshold-level = i] [
      set Identity-group agents with [group = [group] of myself]
      if count Identity-group < 3 [set Identity-group agents]
    ]


    ;; Killing centroids without connected agents
    ask centroids [
      let wom who
      if (not any? agents with [group = wom]) [die]
    ]
    set N_centroids count centroids

    ;; Final coloring and killing of centroids
    if centroid_color? [ask agents [set color (5 + 10 * (group - min [who] of centroids))]]
    if killing_centroids? [ask centroids [die]]
  ]

end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;   G O !   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; Main routine
to go
  ;;;; Preparation part ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Redundant conditions which should be avoided -- if the boundary is drawn as constant, then is completely same whether agents vaguely speak or openly listen,
  ;; seme case is for probability speaking of 100%, then it's same whether individual probability is drawn as constant or uniform, result is still same: all agents has probability 100%.
  if avoid-redundancies? and mode = "vaguely-speak" and boundary-drawn = "constant" [stop]
  if avoid-redundancies? and p-speaking-level = 1 and p-speaking-drawn = "uniform" [stop]

  ;; All preparations of agents, globals, avoiding errors etc. in one sub-routine
  prepare-everything-for-the-step


  ;;;; Main part ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ask agents [
    ;; Firstly, we have to determine satisfaction with the neighborhood via 'get-satisfaction' sub-routine
    set Satisfied? get-satisfaction
    ;; NOTE: Updating of 'Satisfied?' in sub-routine 'prepare-everything-for-the-step' leads to run-time errors.

    ;; Mechanism of own opinion or network change -- decision and rossolution
    ;; In case of dissatisfaction agents changes network
    if not satisfied? [change-of-network]

    ;; Now it's implemented, that agents update opinion with probability regulated by
    ;; slider 'dissatisfied_updates_opinion' or value of variable 'Opponents-ratio'
    ;; in cases they are still not satisfied with their network neighborhood
    ;; TO-DO: Discuss this with the team!
    let value ifelse-value (use_opponents_ratio?) [1 - Opponents-ratio][dissatisfied_updates_opinion]
    if Satisfied? or value > random-float 1 [
      ;if not satisfied? [show "I'm not satisfied!!!" show ticks]
      if model = "HK" [change-opinion-HK]
      ;; Note: Now here is only Hegselmann-Krause algorithm, but in the future we might easily employ other algorithms here!

      ;if not satisfied? [show Opponents-ratio show value]
    ]
  ]

  ;; The main algorithm might produce lonely agents, now we connect them to one/more other speaking agent/s
  connect-loners


  ;;;; Final part ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Recoloring patches, agents, computing how model settled down
  updating-patches-and-globals

  tick

  ;; Finishing condition:
  ;; 1) We reached state, where no turtle changes for RECORD-LENGTH steps, i.e. average of MAIN-RECORD (list of averages of turtles/agents RECORD) is 1 or
  ;; 2) We reached number of steps specified in MAX-TICKS
  recording-situation-and-computing-polarisation
  if (mean main-Record = 1 and network-changes <= 5) or ticks = max-ticks [stop]
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;   SUB/PROCEDURES   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to prepare-everything-for-the-step
  ;; Just checking and avoiding runtime errors part of code
  avoiding-run-time-errors

  ;; Before a round erasing indicator of change
  set network-changes 0

  ;; Update group identities via Louvain
  ifelse use_identity? [
    update-l-distances-weights
    if identity_type = "global" [set-group-identities]
    if identity_type = "individual" [set-individual-group-identities] ; and (ticks mod identity_update = 0)
  ][
    ask agents [set Identity-group agents]
  ]

  ;; speaking, coloring, updating and SATISFACTION!!!
  ask agents [preparing-myself]
end


to preparing-myself
  ;; Updating speaking, color and place
  set speak? speaking
  getColor
  getPlace

  ;; storing previous opinion position as 'Last-opinion'
  set Last-opinion Opinion-position
end

;; for computing thresholds and assigning levels to each agent as per drawing parameters
to compute-identity-thresholds

  ;; computing values
  set id_threshold_set n-values identity_levels [0]
  foreach range identity_levels [ i -> set id_threshold_set replace-item i id_threshold_set precision (i / (identity_levels)) 3]
  ;set id_threshold_set replace-item (identity_levels - 1) id_threshold_set 0.999   ;; since 1 throws an error

  ;; now implementing distribution possibilities
  ;;

  ;; splitting into a major and minor group as 80%-20%. The 80% takes on either the highest (right skew) or lowest values of threshold.
  let eighty_percent_size int N-agents * 0.8
  let major_partition n-of eighty_percent_size agents
  let minor_partition agents with [ not member? self major_partition]
  ifelse (draw_id_threshold = "uniform")[
    ask agents [ set id-threshold-level random identity_levels ]
  ][
    ifelse (identity_levels > 3)[
      ;; if identity levels > 3, only the highest two or lowest two levels will contain 80% of the population
      ask major_partition [
        let dice random 2
        set id-threshold-level dice
        if(draw_id_threshold = "right-skewed") [ set id-threshold-level (id-threshold-level + identity_levels - 2) ]
      ]
      ask minor_partition [
        let dice random (identity_levels - 2)
        set id-threshold-level dice
        if(draw_id_threshold = "left-skewed") [ set id-threshold-level (id-threshold-level + 2) ]
      ]
    ][

      ;; if identity levels < 4, only use one level for major partition

      ask major_partition [
        if(draw_id_threshold = "right-skewed") [ set id-threshold-level identity_levels ]
        if(draw_id_threshold = "left-skewed") [ set id-threshold-level 0 ]
      ]
      ask minor_partition [
        let dice random (identity_levels - 1)
        set id-threshold-level dice
        if(draw_id_threshold = "left-skewed") [ set id-threshold-level (id-threshold-level + 1) ]
      ]
    ]
  ]

end


;; sub-routine computing whether the agent is satisfied with the neighborhood
to-report get-satisfaction
  ;; initialization of agent set 'opponents' containing agents whose opinion positions are key for agent's satisfaction
  let opponents nobody
  ;; 1) updating agent uses only visible link neighbors in small-world network
  let visibles comm-neighbors with [speak?]
  ;; TO-DO: Discuss with Mike, Ashwin and Ashley what scope in Identity group and what in comm-neis.
  ;;        Now I use the IG only for updating opinion in HK, not used in rewiring of comm network, not for satisfaction etc.

  ;; 2) we have different modes for finding supporters:
  ;; 2.1) in mode "I got them!" agent looks outside her boundary (opinion +/- uncertainty),
  ;;      i.e. agent takes as opponents agents that are too far from her opinion
  if mode = "openly-listen" [
    ;; we compute 'lim-dist' -- it is the numerical distance in given opinion space
    let lim-dist (Uncertainty * sqrt(opinions * 4))
    ;; we set as opponents agents with opinion further than 'lim-dist'
    set opponents visibles with [opinion-distance > lim-dist]
  ]

  ;; 2.1) in mode "They got me!" agent looks outside whose boundaries (opinion +/- uncertainty) she is,
  ;;       i.e. she takes as opponents agents that speak with such a low uncertainty that it doesn't match her own opinion
  if mode = "vaguely-speak"  [
    ;; Note: Here is used the 'Uncetainty' value of called agent, agent who might be used for updating,
    ;;       not 'Uncertainty' of calling agent who updates her own opinion.
    set opponents visibles with [opinion-distance > (Uncertainty * sqrt(opinions * 4))]
  ]
  ;; Now we can return the True/False value, whether the agent is satisfied
  ;; (among visible network neighbors from identity group are not too much opponents)
  set Opponents-ratio ifelse-value (count visibles > 0) [precision (count opponents / count visibles) 3][0]  ;; In case no visibles are in the neighborhood, no opponent is there.
  report ifelse-value (count visibles > 0) [Opponents-ratio <= Tolerance][TRUE]  ;; In case no visibles are in the neighborhood, then agent is happy.
end


;; envelope controlling the way, how we change the network
to change-of-network
  if network-change = "link" [rewire-one-link]
  if network-change = "community" [leave-the-neighborhood-join-a-new-one]
  ;; Note: here might be other ways in the future, that's why the 'ifelse' structure is not used here

  ;; We advance the counter of network changes -- now, one just happened.
  set network-changes network-changes + 1

  ;; We check if the changing of network helps
  set Satisfied? get-satisfaction
end


;; subroutine for leaving the neighborhood and joining a new one -- agent is decided to leave, we just process it here
to leave-the-neighborhood-join-a-new-one
  ;; Firstly, we have to count agents neighbors, to determine how many links agent has to create in the main part of the procedure
  let to-visibles my-comms with [[speak?] of other-end]
  let nei-size count to-visibles

  ;; Secondly, we cut off all the links
  ask to-visibles [die]

  ;; Catching possible error with not enough visible agents for creating 'comms'
  let speaking-others other agents with [speak?]
  if (nei-size > count speaking-others) [set nei-size count speaking-others]

  ;; Thirdly, random VS intentional construction of new neighborhood.
  ifelse create-links-randomly? [
    ;; We set new neighborhood randomly or...
    create-comms-with n-of nei-size speaking-others
  ][
    ;; ...creates it out of the closest neighbors.
    create-comms-with min-n-of nei-size speaking-others [opinion-distance]
  ]

  ;; P.S. Just hiding links for better speed of code -- when we change/cut a link, all links become visible and that slows down the simulation.
  ask comms [set hidden? TRUE]
end

;; TO-DO: agents should cut-off only neighbors that they previously heard speak,
;;        we probably should create their memory whom they heard speak and onlythose agents might cut-off.
;;
;; Note:  Now I implement it in modest variant: agent cuts off the most annoying/random presently speaking agent --
;;        there must be at least one, since agents are satisfied by the rule with the empty neighborhood and
;;        they update neighborhood only in case of dissatisfaction.

;; subroutine for changing one link
to rewire-one-link
  ;; Firstly, we find speaking agents who are both: in Identity group and connected by communication link.
  let visibles comm-neighbors with [speak?]

  ;; Cutting-off of the link itself:
  let a-visible one-of visibles
  let annoyer max-one-of visibles [opinion-distance]
  ask one-of my-comms with [other-end = ifelse-value (cut-links-randomly?) [a-visible][annoyer]] [die]

  ;; Secondly, we choose for the agent a new speaking partner with the random/most-close opinion
  ;let cha count my-comms
  let potentials other agents with [speak? and not comm-neighbor? myself]
  ;if 128 != (count potentials + cha) [show "!!!!!!!!!!!!!!!!!!!MISTAKE!!!!!!!!!!!!!!!!!!!! wrong potentials/my-comms" show ticks]
  ;if count potentials = 0 [show "!!!!!!!!!!!!!!!!!!!MISTAKE!!!!!!!!!!!!!!!!!!!! zero potentials" show ticks]
  let a-partner one-of potentials
  let the-partner min-one-of potentials [opinion-distance]
  create-comm-with ifelse-value (create-links-randomly?) [a-partner] [the-partner]
  ;if not (count my-comms > cha) [show "!!!!!!!!!!!!!!!!!!!MISTAKE!!!!!!!!!!!!!!!!!!!! no link created!" show ticks]

  ;; P.S. Just hiding links for better speed of code -- when we change/cut a link, all links become visible and that slows down the simulation.
  ask comms [set hidden? TRUE]
end


;; Procedure for connecting agents with not enough comm-neigbours
to connect-loners
  ;; We check whether each agent has enough neighbors
  ask agents with [count comm-neighbors < min-comm-neis] [
    ;; !!!WARNING!!!: If there are many agents classified as 'loners',
    ;; then early running agents might connect some later going agents,
    ;; and then might happen that when later running agent runs this algorithm,
    ;; she has enough connections (because she was connected by earlier running agents).
    ;; SOLUTION: Check before creating new link whetwer agent is still demanding these new links.

    ;; Defining needed agentset and variables:
    let potentials other agents with [speak?]  ; We set 'potentials' to all other speaking agents and then...
    let p count potentials  ; 'p' stands for potentials

    ;; Catch of potential BUG via 'if' structure --
    ;;   a) if there is no-one speaking, then the lone agent has to wait until the next round.
    ;;   b) if agent was demanding new links, but was served by previous demanders, then needs no new link
    if p > 0 and count comm-neighbors < min-comm-neis [
        let n min-comm-neis - count comm-neighbors  ; 'n' stands for needed
        let ap ifelse-value(p >= n)[n][p]  ; 'ap' stands for asked potentials
        create-comms-with ifelse-value (create-links-randomly?) [n-of ap potentials][min-n-of ap potentials [opinion-distance]]  ;... it depends on scenario: we choose randomly or with the closest opinion
    ]
  ]

  ;; P.S. Just hiding links for better speed of code -- when we change/cut a link, all links become visible and that slows down the simulation.
  ask comms [set hidden? TRUE]
end


;; sub-routine for updating opinion position of turtle according the Hegselmann-Krause (2002) model
to change-opinion-HK
  ;; initialization of agent set 'influentials' containing agents whose opinion positions uses updating agent
  let influentials no-turtles
  ;; 1) updating agent uses only visible link neighbors in small-world network who are also members of her Identity group
  let visibles (turtle-set filter [vis -> member? vis Identity-Group] sort comm-neighbors with [speak?])

  ;; 2) we have different modes for finding influentials:
  ;; 2.1) in mode "openly-listen" agent looks inside his boundary (opinion +/- uncertainty),
  ;;      i.e. agent takes opinions not that much far from her opinion
  if mode = "openly-listen" [
    ;; we compute 'lim-dist' -- it is the numerical distance in given opinion space
    let lim-dist (Uncertainty * sqrt(opinions * 4))
    ;; we set as influentials agents with opinion not further than 'lim-dist'
    set influentials visibles with [opinion-distance <= lim-dist]
  ]

  ;; 2.1) in mode "vaguely-speak" agent looks inside whose boundaries (opinion +/- uncertainty)
  ;;       she is, i.e. agents takes opinions spoken with such a big uncertainty that it matches her own opinion
  if mode = "vaguely-speak"  [
    ;; Note: Here is used the 'Uncetainty' value of called agent, agent who might be used for updating,
    ;;       not 'Uncertainty' of calling agent who updates her opinion.
    set influentials visibles with [opinion-distance <= (Uncertainty * sqrt(opinions * 4))]
  ]

  ;; 3) we also add the updating agent into 'influentials'
  set influentials (turtle-set self influentials)

  ;; we check whether there is someone else then calling/updating agent in the agent set 'influentials'
  if count influentials > 1 [

    ;; here we draw a list of dimensions which we will update:
    ;; by 'range opinions' we generate list of integers from '0' to 'opinions - 1',
    ;; by 'n-of updating' we randomly take 'updating' number of integers from this list
    ;; by 'shuffle' we randomize order of the resulting list
    let op-list shuffle n-of updating range opinions

    ;; we initialize counter 'step'
    let step 0

    ;; we go through the while-loop 'updating' times:
    while [step < updating] [
      ;; we initialize/set index of updated opinion dimension according the items on the 'op-list',
      ;; note: since we use while-loop, we go through each item of the 'op-list', step by step, iteration by iteration.
      let i (item step op-list)

      ;; then we update dimension of index 'i' drawn from the 'op-list' in the previous line:
      ;; 1) we compute average position in given dimension of the calling/updating agent and all other agents from agent set 'influentials'
      ;;    by the command '(mean [item i opinion-position] of influentials)', and
      ;; 2) the new value of opinion 'val' is not directly average, but it is weighted by the 'Conformity' (individual trait),
      ;;    the closer 'Conformity' to 1, the closer agent jumps into the mean of others, the closer to 0, the less agent moves.
      ;; 3) we set value as new opinion position by command 'set opinion-position replace-item i opinion-position X' where 'X' is the mean opinion (ad 1, see line above)

      ;; ad 1: averge position computation
      let val  precision (mean [item i opinion-position] of influentials) 3 ;; NOTE: H-K model really assumes that agent adopts immediatelly the 'consesual' position
      ;; ad 2: updating/weighting 'val' by 'Conformity' and own opinion
      let my item i opinion-position
      set val my + ((val - my) * Conformity)
      ;; ad 3: assigning the value 'val'
      set opinion-position replace-item i opinion-position val

      ;; advancement of counter 'step'
      set step step + 1
    ]
  ]
end


;; sub-routine for computing opinion distance of two comparing agents
to-report opinion-distance
  ;; we store in temporary variable the opinion of the called and compared agent
  let my opinion-position

  ;; we store in temporary variable the opinion of the calling and comparing agent
  let her [opinion-position] of myself

  ;; we initialize counter of step of comparison -- we will compare as many times as we have dimensions
  let step 0

  ;; we initialize container where we will store squared distance in each dimension
  let dist 0

  ;; while loop going through each dimension, computiong distance in each dimension, squarring it and adding in the container
  while [step < opinions] [
    ;; computiong distance in each dimension, squarring it and adding in the container
    set dist dist + (item step my - item step her) ^ 2

    ;; advancing 'step' counter by 1
    set step step + 1
  ]

  ;; computing square-root of the container 'dist' -- computing Euclidean distance -- and setting it as 'dist'
  set dist sqrt dist

  ;; reporting Euclidean distance
  report dist
end


;; sub-routine for computing opinion distance of two comparing opinion positions  -- relative distance weighted as 1 for minimal distance and 0 for the maximal one
to-report opinion-distance2 [my her]
  ;; we initialize counter of step of comparison -- we will compare as many times as we have dimensions
  let step 0

  ;; we initialize container where we will store squared distance in each dimension
  let dist 0

  ;; while loop going through each dimension, computiong distance in each dimension, squarring it and adding in the container
  while [step < opinions] [
    ;; computiong distance in each dimension, squarring it and adding in the container
    set dist dist + (item step my - item step her) ^ 2

    ;; advancing 'step' counter by 1
    set step step + 1
  ]

  ;; computing square-root of the container 'dist' -- computing Euclidean distance -- and setting it as 'dist'
  set dist sqrt dist

  ;; Turning 'dist' into 'weight'
  let weight (sqrt(4 * opinions) - dist) / sqrt(4 * opinions)

  ;; reporting weight of distance
  report precision weight 10
end


;; sub-routine for computing opinion distance of two comparing opinion positions  -- absolute distance without weighting
to-report opinion-distance3 [my her]
  ;; we initialize counter of step of comparison -- we will compare as many times as we have dimensions
  let step 0

  ;; we initialize container where we will store squared distance in each dimension
  let dist 0

  ;; while loop going through each dimension, computiong distance in each dimension, squarring it and adding in the container
  while [step < opinions] [
    ;; computiong distance in each dimension, squarring it and adding in the container
    set dist dist + (item step my - item step her) ^ 2

    ;; advancing 'step' counter by 1
    set step step + 1
  ]

  ;; computing square-root of the container 'dist' -- computing Euclidean distance -- and setting it as 'dist'
  set dist sqrt dist

  ;; reporting weight of distance
  report precision dist 10
end


to recording-situation-and-computing-polarisation
  ;; Finishing condition:
  ;; 1) We reached state, where no turtle changes for RECORD-LENGTH steps, i.e. average of MAIN-RECORD (list of averages of turtles/agents RECORD) is 1 or
  ;; 2) We reached number of steps specified in MAX-TICKS
  if ((mean main-Record = 1 and network-changes <= 5) or ticks = max-ticks) [compute-polarisation-repeatedly]
  if ((mean main-Record = 1 and network-changes <= 5) or ticks = max-ticks) [compute-final-macro-state-of-simulation]

  ;; Recording and computing polarisation on the fly...
  ;if (ticks / polarisation-each-n-steps) = floor (ticks / polarisation-each-n-steps) [compute-polarisation-repeatedly]
  ;if (ticks / record-each-n-steps) = floor(ticks / record-each-n-steps) [compute-final-macro-state-of-simulation]
end


;; Updating patches and global variables
to updating-patches-and-globals
  ;; Patches update color according the number of turtles on it.
  ask patches [set pcolor patch-color]

  ;; We have to check here the change of opinions, resp. how many agents changed,
  ;; and record it for each agent and also for the whole simulation
  ;; Turtles update their record of changes:
  ask agents [
    ;; we take 1 if opinion is same, we take 0 if opinion changes, then
    ;; we put 1/0 on the start of the list Record, but we omit the last item from Record
    set Record fput ifelse-value (Last-opinion = Opinion-position) [1][0] but-last Record
  ]
  ;; Then we might update it for the whole:
  set main-Record fput precision (mean [mean Record] of agents) 3 but-last main-Record

  ;; Coloring agents according identity group
  if centroid_color? [ask agents [set color (5 + 10 * (group - min [group] of agents))]]
end


;; Procedure reporting ESBG/Ashwin's polarisation
to-report Ash-polarisation
  ;; Preparation
  create-centroids 2 [set shape "square" set Opinion-position n-values opinions [precision (1 - random-float 2) 8]]
  let cent1 max [who] of centroids  ;; Storing 'who' of two new centroids
  let cent0 cent1 - 1
  ask agents [set group (cent0 + (who mod 2))]  ;; Random assignment of agents to the groups
  updating-centroids-opinion-position (cent0) (cent1)  ;; Initial update

  ;; Iterating until centroids are stable
  while [Centroids_change < sum [opinion-distance3 Opinion-position Last-opinion] of centroids with [who >= cent0]][
    update-agents-opinion-group (cent0) (cent1)
    updating-centroids-opinion-position (cent0) (cent1)
  ]

  ;; Computing polarisation -- cutting-out agents too distant from centroids
  ask agents [set distance_to_centroid [opinion-distance] of centroid group]
  let a0 agents with [group = cent0]
  set a0 min-n-of (count a0 - ESBG_furthest_out) a0 [opinion-distance3 ([opinion-position] of self) ([opinion-position] of centroid cent0)]
  let a1 agents with [group = cent1]
  set a1 min-n-of (count a1 - ESBG_furthest_out) a1 [opinion-distance3 ([opinion-position] of self) ([opinion-position] of centroid cent1)]

  ;; Updating centroids and agents opinion position (without furthest agents)
  ask centroid cent0 [foreach range opinions [o -> set opinion-position replace-item o opinion-position precision (mean [item o opinion-position] of a0) 8] getPlace]
  ask a0 [set distance_to_centroid [opinion-distance] of centroid cent0]
  ask centroid cent1 [foreach range opinions [o -> set opinion-position replace-item o opinion-position precision (mean [item o opinion-position] of a1) 8] getPlace]
  ask a1 [set distance_to_centroid [opinion-distance] of centroid cent1]

  ;; Preparing final distances and diversity
  let cent-dist opinion-distance3 ([opinion-position] of centroid cent0) ([opinion-position] of centroid cent1)
  let div0 (mean [distance_to_centroid] of a0)
  let div1 (mean [distance_to_centroid] of a1)

  ;; Cleaning and reporting
  ask centroids with [who >= cent0] [die]
  report (cent-dist / (1 + div0 + div1)) / sqrt(opinions * 4)
end


to update-agents-opinion-group [cent0 cent1]
  ;; Checking the assignment -- is the assigned centroid the nearest? If not, reassign!
  ask agents [set group group - ([who] of min-one-of centroids with [who >= cent0] [opinion-distance])]  ; set color 15 + group * 10]
  let wrongly-at-grp0 turtle-set agents with [group = -1]  ;; they are in 0, but should be in 1: 0 - 1 = -1
  let wrongly-at-grp1 turtle-set agents with [group = 1]  ;; they are in 1, but should be in 0: 1 - 0 = 1
  ifelse count wrongly-at-grp0 = count wrongly-at-grp1 [
    ask agents [set group [who] of min-one-of centroids with [who >= cent0] [opinion-distance]]
  ][
    let peleton agents with [group = 0]
    ifelse count wrongly-at-grp0 < count wrongly-at-grp1 [
      set peleton (turtle-set peleton wrongly-at-grp0 max-n-of (count wrongly-at-grp0) wrongly-at-grp1 [opinion-distance3 ([opinion-position] of self) ([opinion-position] of centroid cent0)]) ;; all agents assigned correctly + smaller group of wrong + from bigger group 'n of size of smaller group'
      let stayed agents with [not member? self peleton]
      ask peleton [set group [who] of min-one-of centroids with [who >= cent0] [opinion-distance] ;set color 15 + group * 10
      ]
      ask stayed [set group cent1 ;set color 15 + group * 10
      ]
     ][
      set peleton (turtle-set peleton wrongly-at-grp1 max-n-of (count wrongly-at-grp1) wrongly-at-grp0 [opinion-distance3 ([opinion-position] of self) ([opinion-position] of centroid cent1)]) ;; all agents assigned correctly + smaller group of wrong + from bigger group 'n of size of smaller group'
      let stayed agents with [not member? self peleton]
      ask peleton [set group [who] of min-one-of centroids with [who >= cent0] [opinion-distance] ;set color 15 + group * 10
      ]
      ask stayed [set group cent0 ;set color 15 + group * 10
      ]
    ]
  ]
end


to updating-centroids-opinion-position [cent0 cent1]
  ;; Storing opinion as Last-opinion
  ask centroids with [who >= cent0] [set Last-opinion Opinion-position]

  ;; Computing groups mean 'opinion-position'
  set positions_clusters [] ;; List with all positions of both 2 groups
  foreach range 2 [c ->
    let one-position []  ;; List for one position of one cluster
    foreach range opinions [o -> set one-position lput precision (mean [item o opinion-position] of agents with [group = cent0 + c]) 8 one-position]
    set positions_clusters lput one-position positions_clusters
  ]

  ;; Setting opinions of centroids
  ask centroid cent0 [set opinion-position item 0 positions_clusters getPlace]
  ask centroid cent1 [set opinion-position item 1 positions_clusters getPlace]
end


;; Envelope for combined update of weights of all links
to update-links-weights
  update-comms-weights
  update-l-distances-weights
end


;; Sub-routine for updating Communication links' weights,
;; according opinion distance of both their ends
to update-comms-weights
  ;; We use function 'opinion-distance2', which needs two opinion positions as input and
  ;; receives their distance as output, but this distance is converted to weight:
  ;; weight = 1 means that both positions are same, weight = 0 means that their distance is maximal,
  ;; i.e. both positions are in oposit corners of respective N-dimensional space.
  ask comms [set op-weight opinion-distance2 ([opinion-position] of end1) ([opinion-position] of end2)]
  ask comms [set hidden? TRUE]
end


;; Sub-routine for updating Distances links' weights,
;; according opinion distance of both their ends
to update-l-distances-weights
  ;; We use function 'opinion-distance2', which needs two opinion positions as input and
  ;; receives their distance as output, but this distance is converted to weight:
  ;; weight = 1 means that both positions are same, weight = 0 means that their distance is maximal,
  ;; i.e. both positions are in oposit corners of respective N-dimensional space.
  ask l-distances [set l-weight opinion-distance2 ([opinion-position] of end1) ([opinion-position] of end2)]
  ask l-distances [set hidden? TRUE]
end


;; We compute polarisation several times and then set it for the average
to compute-polarisation-repeatedly
  ;; Initialization of temporal variables
  let r 0
  let p []
  let np []
  let up []
  let unp []
  let ap []
  update-links-weights

  ;; Repeating cycle
  while [r < polar_repeats] [
    compute-polarisation
    set p lput polarisation p
    set np lput normalized_polarisation np
    set up lput unweighted_polarisation up
    set unp lput unweighted_normalized_polarisation unp
    set ap lput Ash-polarisation ap
    set r r + 1
  ]

  ;; Setting variables back
  set polarisation precision (mean p) 3
  set normalized_polarisation precision (mean np) 3
  set unweighted_polarisation precision (mean up) 3
  set unweighted_normalized_polarisation precision (mean unp) 3
  set ESBG_polarisation precision (mean ap) 3
end

;; NOTE: Now I am iplementing it for N = 2 centroids, but I prepare code for easy generalisation for N > 2.
to compute-polarisation
  ;; Preparation
  ask centroids [die]

  ;; Detection of clusters via Louvain
  let selected-agents agents with [2 < count my-l-distances with [l-weight >= d_threshold]]  ;; Note: We take into account only not loosely connected agents
  nw:set-context selected-agents l-distances with [l-weight >= d_threshold]
  let communities nw:louvain-communities
  set N_centroids length communities

  ;; Computing clusters' mean 'opinion-position'
  let positions-clusters [] ;; List with all positions of all clusters
  foreach communities [c ->
    let one []  ;; List for one positio nof one cluster
    foreach range opinions [o -> set one lput precision (mean [item o opinion-position] of c) 8 one]
    set positions-clusters lput one positions-clusters
  ]

  ;; Preparation of centroids -- feedeing them with communities
  create-centroids N_centroids [
    set heading (who - min [who] of centroids)
    set Opinion-position item heading positions-clusters  ;; We set opinions, we try to do it smoothly...
    set shape "circle"
    set size 1.5
    set color 5 + (who - min [who] of centroids) * 10
    getPlace
  ]

  ;; Assignment of agents to groups
  ask selected-agents [set group [who] of min-one-of centroids [opinion-distance]]

  ;; Computation of centroids possitions
  compute-centroids-positions (selected-agents)

  ;; Iterating cycle -- looking for good match of centroids
  while [sum [opinion-distance3 (Last-opinion) (Opinion-position)] of centroids > Centroids_change] [

    ;; turtles compute whether they are in right cluster and
    ask selected-agents [set group [who] of min-one-of centroids [opinion-distance]]

    ;; Computation of centroids possitions
    compute-centroids-positions (selected-agents)
  ]

  ;; Killing centroids without connected agents
  ask centroids [
    let wom who
    if (not any? agents with [group = wom]) [die]
  ]
  set N_centroids count centroids

  ;; Catching the run-time error:
  ;; If there is just one component since all agents has same opinion, then the polarisation algorithm does produce error --
  ;; because of computing mean of empty list of distences: this list is empty since the only one existing centroid can't
  ;; compute distance to itself via double 'while' structuresince 'ai' and 'aj' lists are empty.
  ;; In this case it is obvious that polarisation is 0, so we set 'polarisation' and 'normalized_polarisation' to 0 directly via 'ifelse' structure
  ifelse (count centroids < 2) [
    set polarisation 0
    set normalized_polarisation 0
    set unweighted_polarisation 0
    set unweighted_normalized_polarisation 0
  ][
    ;; Computing polarization -- preparation of lists and agents
    let distances []
    let diversity []
    let unweighted_distances []
    let unweighted_diversity []
    let whos sort [who] of centroids  ;; List of 'who' of all centroids
    ask selected-agents [set distance_to_centroid [opinion-distance] of centroid group]  ;; Each agent computes her distance to her centroid and stores it as 'distance_to_centroid'.

    ;; Computing polarization -- distances of all centroids
    let ai but-last whos  ;; List of all 'i' -- 'who' initializing distances computation
    let aj but-first whos  ;; List of all 'j' -- 'who' of other end of distances computation
    foreach ai [i ->
      foreach aj [j ->
        ;; Each distance is weighted by fraction of both centroid groups
        ;; via formula '(N_centroids ^ 2) * (count agents with [group = i] / count agents) * (count agents with [group = j] / count agents)'
        let weight (N_centroids ^ 2) * (count selected-agents with [group = i] / count selected-agents) * (count selected-agents with [group = j] / count selected-agents)
        let cent-dist opinion-distance3 ([opinion-position] of centroid i) ([opinion-position] of centroid j)
        set distances lput (weight * cent-dist) distances
        set unweighted_distances lput cent-dist unweighted_distances
      ]
      set aj but-first aj
    ]

    ;; Computing polarization -- diversity in groups
    foreach sort [who] of centroids [wg ->
      let weight (count selected-agents with [group = wg] / count selected-agents)
      let cent-div (mean [distance_to_centroid] of selected-agents with [group = wg])
      set diversity lput (weight * cent-div) diversity
      set unweighted_diversity lput cent-div unweighted_diversity
    ]

    ;; Computing polarization -- polarization indexes
    ;; Note: Now it is computed to receive same result as for n=2 groups,
    ;;       but might be needed to change it later to get better results for n>2 group,
    ;;       but now it works fine without runtime errors with all numbers of groups.
    set polarisation (mean distances) / (1 + 2 * mean diversity)  ;; Raw polarization computed as distance divided by heterogeinity in the groups.
    set normalized_polarisation precision (polarisation / (2 * sqrt(opinions))) 3
    set polarisation precision polarisation 3
    set unweighted_polarisation (mean unweighted_distances) / (1 + 2 * mean unweighted_diversity)  ;; Raw unweighted polarization computed as unweighted distance divided by unweighted heterogeinity in the groups.
    set unweighted_normalized_polarisation precision (unweighted_polarisation / (2 * sqrt(opinions))) 3
    set unweighted_polarisation precision unweighted_polarisation 3
  ]

  ;; Final coloring and killing of centroids
  if centroid_color? [ask selected-agents [set color (5 + 10 * (group - min [who] of centroids))]]
  if killing_centroids? [ask centroids [die]]
end


;; Sub-routine of polarization routine
to compute-centroids-positions [sel-agents]
  ;; Preparation
  ask centroids [set Last-opinion Opinion-position]

  ;; Computation of centroids possitions
  let grp min [who] of centroids
  while [grp <= (max [who] of centroids)] [
    ask centroid grp [
      ifelse (not any? agents with [group = grp]) [
        set Opinion-position Last-opinion
      ][
        let dim 0
        while [dim < opinions] [
          set Opinion-position replace-item dim Opinion-position mean [item dim Opinion-position] of sel-agents with [group = grp]
          set dim dim + 1
        ]
      ]
    ]
    set grp grp + 1
  ]
  ask centroids [
    getPlace
  ]
end


;; reporter function for translating a list into one string of values divided by commas
to-report list-to-string [LtS]
  ;; Initializing empty string and counter
  let str ""
  let i 0

  ;; Now we go through the list item by item and add them into string
  while [i < length LtS][
    set str (word str item i LtS)
    set i i + 1
    if (i < length LtS) [set str (word str ", ")]
  ]

  report str
end


;; Sub-routine for assigning value of tolerance
to-report get-group-threshold
  ;; We have to initialize empty temporary variable
  let gtValue 0

  ;; Then we draw the value according the chosen method
  if threshold_drawn = "constant" [set gtValue id_threshold + random-float 0]  ;; NOTE! 'random-float 0' is here for consuming one pseudorandom number to cunsume same number of pseudorandom numbers as "uniform
  if threshold_drawn = "uniform" [set gtValue ifelse-value (id_threshold < 0.5)
                                                           [precision (random-float (2 * id_threshold)) 3]
                                                           [precision (1 - (random-float (2 * (1 - id_threshold)))) 3]]
  report gtValue
end


;; Sub-routine for assigning value of tolerance
to-report get-conformity
  ;; We have to initialize empty temporary variable
  let cValue 0

  ;; Then we draw the value according the chosen method
  if conformity-drawn = "constant" [set cValue conformity-level + random-float 0]  ;; NOTE! 'random-float 0' is here for consuming one pseudorandom number to cunsume same number of pseudorandom numbers as "uniform
  if conformity-drawn = "uniform" [set cValue ifelse-value (conformity-level < 0.5)
                                                           [precision (random-float (2 * conformity-level)) 3]
                                                           [precision (1 - (random-float (2 * (1 - conformity-level)))) 3]]
  report cValue
end

;; Sub-routine for assigning value of tolerance
to-report get-tolerance
  ;; We have to initialize empty temporary variable
  let tValue 0

  ;; Then we draw the value according the chosen method
  if tolerance-drawn = "constant" [set tValue tolerance-level + random-float 0]  ;; NOTE! 'random-float 0' is here for consuming one pseudorandom number to cunsume same number of pseudorandom numbers as "uniform
  if tolerance-drawn = "uniform" [set tValue ifelse-value (tolerance-level < 0.5)
                                                           [precision (random-float (2 * tolerance-level)) 3]
                                                           [precision (1 - (random-float (2 * (1 - tolerance-level)))) 3]]
  report tValue
end


;; Sub-routine for assigning value of p-speaking
to-report get-speaking
  ;; We have to initialize empty temporary variable
  let pValue 0

  ;; Then we draw the value according the chosen method
  if p-speaking-drawn = "constant" [set pValue p-speaking-level + random-float 0]  ;; NOTE! 'random-float 0' is here for consuming one pseudorandom number to cunsume same number of pseudorandom numbers as "uniform
  if p-speaking-drawn = "uniform" [set pValue ifelse-value (p-speaking-level < 0.5)
                                                           [precision (random-float (2 * p-speaking-level)) 3]
                                                           [precision (1 - (random-float (2 * (1 - p-speaking-level)))) 3]]
  if p-speaking-drawn = "function" [set pValue (precision(sqrt (sum (map [ x -> x * x ] opinion-position)) / sqrt opinions) 3) + random-float 0]  ;; NOTE! 'random-float 0' is here for consuming one pseudorandom number to cunsume same number of pseudorandom numbers as "uniform

  ;; Report result back
  report pValue
end


;; sub-routine for assigning value of uncertainty to the agent
to-report get-uncertainty
  ;; We have to initialize empty temporary variable
  let uValue 0

  ;; Then we draw the value according the chosen method
  if boundary-drawn = "constant" [set uValue boundary + random-float 0]  ;; NOTE! 'random-float 0' is here for consuming one pseudorandom number to cunsume same number of pseudorandom numbers as "uniform"
  if boundary-drawn = "uniform" [set uValue precision (random-float (2 * boundary)) 3]

  ;; reporting value back for assigning
  report uValue
end


;; sub-routine for graphical representation -- it takes two opinion dimension and gives the agent on XY coordinates accordingly
to getPlace
  ;; check whether our cosen dimension is not bigger than maximum of dimensions in the simulation
  if X-opinion > opinions [set X-opinion 1]
  if Y-opinion > opinions [set Y-opinion 1]

  ;; then we rotate the agent towards the future place
  facexy ((item (X-opinion - 1) opinion-position) * max-pxcor) ((item (Y-opinion - 1) opinion-position) * max-pycor)

  ;; lastly we move agent on the place given by opinion dimensions chosen for X and Y coordinates
  set xcor (item (X-opinion - 1) opinion-position) * max-pxcor
  set ycor (item (Y-opinion - 1) opinion-position) * max-pycor
end


;; sub routine for coloring agents according their average opinion across all dimensions --
;; useful for distinguishing agents with same displayed coordinates, but differing in other opinion dimensions,
;; then we see at one place agents with different colors.
to getColor
  ;; speaking agents are colored from very dark red (average -1) through red (average 0) to very light red (average +1)
  ifelse speak? [
    set color 15 + 4 * mean(opinion-position)
    set size (max-pxcor / 10)
  ]
  ;; silent agent are white and of zero size, to just show the speaking one -- later we might parametrize this if we want...
  [
    set color white
    set size 0
  ]
end


;; Sub routine for dissolving whether agent speaks at the given round/step or not
to-report speaking
  ;; For the case of function we have to update P-speaking value
  if p-speaking-drawn = "function" [set P-speaking precision(sqrt (sum (map [ x -> x * x ] opinion-position)) / sqrt opinions) 3]

  report P-speaking > random-float 1
end


;; sub-routine for visual purposes -- colors empty patches white, patches with some agents light green, with many agents dark green, with all agents black
to-report patch-color
  report 59.9 - (9.8 * (ln(1 + count turtles-here) / ln(N-agents)))
end


;; Sub routine just for catching run-time errors
to avoiding-run-time-errors
  ;; Check whether we set properly parameter 'updating' --
  ;; if we want update more dimensions than exists in simulation, then we set 'updating' to max of dimensions, i.e. 'opinions'
  if updating > opinions [set updating opinions]
end


;; Sub-routine which opens/creates *.csv file and stores there states of all turtles
to record-state-of-simulation
  ;; setting working directory
  ;set-current-directory directory

  ;; seting 'file-name'
  let file-name (word "Sims/Nodes01_" file-name-core "_" ticks ".csv")

  ;;;; File creation and opening: NODES
  ;; If file exists at the start we delete it to start with clean file
  if file-exists? file-name [file-delete file-name]
  file-open file-name ;; This opens existing file (at the end) or creates file if doesn't exist (at the begining)

    ;; Constructing list for the first row with variable names:
    let row (list "ID" "Uncertainty" "pSpeaking" "Speaks")
    foreach range Opinions [i -> set row lput (word "Opinion" (i + 1)) row]

    ;; Writing the variable names in the first row at the start
    file-print list-to-string (row)

  ;; For writing states itself we firstly need to create list of sorted turtles 'srt'
  let srt sort agents
  ;; Then we iterate over the list 'srt':
  foreach srt [t -> ask t [  ;; every turtle in the list...
    set row (list (word "Nei" who)  Uncertainty P-Speaking (ifelse-value (speak?)[1][0]))  ;; stores in list ROW its ID, Uncertainty, P-Speaking and whether speaks...
    foreach Opinion-position [op -> set row lput (precision(op) 3) row]  ;; Opinions ...
    file-print list-to-string (row)
    file-flush  ;; for larger simulations with many agents it will be safer flush the file after each row
  ]]

  ;; Finally, we close the file
  file-close
  file-close


  ;;;; File creation and opening: LINKS
  ;; seting 'file-name' for links.
  set file-name (word "Sims/Links01_" file-name-core "_" ticks ".csv")

  ;;;; File creation and opening
  ;; If file exists at the start we delete it to start with clean file
  if file-exists? file-name [file-delete file-name]
  file-open file-name ;; This opens existing file (at the end) or creates file if doesn't exist (at the begining)

    ;; Constructing list for the first row with variable names:
    set row (list "ID1" "ID2" "Communication" "Distance")

    ;; Writing the variable names in the first row at the start
    file-print list-to-string (row)

  ;; We need to prepare counters and other auxilliary varibles for doubled cycle:
  let i 0
  let j 1
  let iMax (count agents - 2)
  let jMax (count agents - 1)
  let mine []
  let her []

  ;; Now double while cycle!
  while [i <= iMax] [  ;; Iterating over all turtles except the last one
    set j i + 1
    while [j <= jMax][  ;; Second cycle iterates over all turtles with index higher than 'i'
      set mine ([opinion-position] of turtle i) ;; First opinion for measuring distance...
      set her ([opinion-position] of turtle j)  ;; Second opinion...
      set row (list i j (ifelse-value (is-comm? comm i j) [1][0]) opinion-distance2 (mine) (her))  ;; Construction of the 'row'
      file-print list-to-string (row)  ;; Writing the row...
      file-flush  ;; for larger simulations with many agents it will be safer flush the file after each row
      set j j + 1 ;; Updating counter of second cycle
    ]
    set i i + 1  ;; Updating counter of the first cycle
  ]

  ;; Finally, we close the file
  file-close
  file-close
end


;; Subroutine for computing aggregate/macro state of simulation
to compute-initial-macro-state-of-simulation
    ;; Network
    nw:set-context turtles comms ;; Setting context for the network measures
    set betweenness_start mean [nw:betweenness-centrality] of turtles
    set eigenvector_start mean [nw:eigenvector-centrality] of turtles
    set clustering_start mean [nw:clustering-coefficient] of turtles
    ;modularity_start
    set mean_path_start nw:mean-path-length

    ;; Polarisation
    compute-polarisation-repeatedly ;; Firstly we have to compute polarization
    set normalized_polarization_start normalized_polarisation
    set ESBSG_polarization_start ESBG_polarisation

    ;; Opinion
    set mean_op1_start mean [item 0 opinion-position] of turtles
    set mean_op2_start mean [item 1 opinion-position] of turtles
    set mean_op3_start mean [item 2 opinion-position] of turtles
    set mean_op4_start mean [item 3 opinion-position] of turtles
    set sd_op1_start standard-deviation [item 0 opinion-position] of turtles
    set sd_op2_start standard-deviation [item 1 opinion-position] of turtles
    set sd_op3_start standard-deviation [item 2 opinion-position] of turtles
    set sd_op4_start standard-deviation [item 3 opinion-position] of turtles

    ;; For last measures we need to sort turtles regarding their opinion
    let lops1 sort-by < [item 0 opinion-position] of turtles
    let lops2 sort-by < [item 1 opinion-position] of turtles
    let lops3 sort-by < [item 2 opinion-position] of turtles
    let lops4 sort-by < [item 3 opinion-position] of turtles
    set median_op1_start item 64 lops1
    set median_op2_start item 64 lops2
    set median_op3_start item 64 lops3
    set median_op4_start item 64 lops4
    set lower_op1_start mean (list item 31 lops1 item 32 lops1)
    set lower_op2_start mean (list item 31 lops2 item 32 lops2)
    set lower_op3_start mean (list item 31 lops1 item 32 lops1)
    set lower_op4_start mean (list item 31 lops2 item 32 lops2)
    set upper_op1_start mean (list item 95 lops1 item 96 lops1)
    set upper_op2_start mean (list item 95 lops2 item 96 lops2)
    set upper_op3_start mean (list item 95 lops1 item 96 lops1)
    set upper_op4_start mean (list item 95 lops2 item 96 lops2)
end
to compute-final-macro-state-of-simulation
    ;; Network
    nw:set-context turtles comms ;; Setting context for the network measures
    set betweenness_final mean [nw:betweenness-centrality] of turtles
    set eigenvector_final mean [nw:eigenvector-centrality] of turtles
    set clustering_final mean [nw:clustering-coefficient] of turtles
    ;modularity_final
    set mean_path_final nw:mean-path-length

    ;; Polarisation
    compute-polarisation-repeatedly ;; Firstly we have to compute polarization
    set normalized_polarization_final normalized_polarisation
    set ESBSG_polarization_final ESBG_polarisation

    ;; Opinion
    set mean_op1_final mean [item 0 opinion-position] of turtles
    set mean_op2_final mean [item 1 opinion-position] of turtles
    set mean_op3_final mean [item 2 opinion-position] of turtles
    set mean_op4_final mean [item 3 opinion-position] of turtles
    set sd_op1_final standard-deviation [item 0 opinion-position] of turtles
    set sd_op2_final standard-deviation [item 1 opinion-position] of turtles
    set sd_op3_final standard-deviation [item 2 opinion-position] of turtles
    set sd_op4_final standard-deviation [item 3 opinion-position] of turtles

    ;; For last measures we need to sort turtles regarding their opinion
    let lops1 sort-by < [item 0 opinion-position] of turtles
    let lops2 sort-by < [item 1 opinion-position] of turtles
    let lops3 sort-by < [item 2 opinion-position] of turtles
    let lops4 sort-by < [item 3 opinion-position] of turtles
    set median_op1_final item 64 lops1
    set median_op2_final item 64 lops2
    set median_op3_final item 64 lops1
    set median_op4_final item 64 lops2
    set lower_op1_final mean (list item 31 lops1 item 32 lops1)
    set lower_op2_final mean (list item 31 lops2 item 32 lops2)
    set lower_op3_final mean (list item 31 lops1 item 32 lops1)
    set lower_op4_final mean (list item 31 lops2 item 32 lops2)
    set upper_op1_final mean (list item 95 lops1 item 96 lops1)
    set upper_op2_final mean (list item 95 lops2 item 96 lops2)
    set upper_op3_final mean (list item 95 lops1 item 96 lops1)
    set upper_op4_final mean (list item 95 lops2 item 96 lops2)
end
@#$#@#$#@
GRAPHICS-WINDOW
219
10
667
459
-1
-1
40.0
1
10
1
1
1
0
0
0
1
-5
5
-5
5
1
1
1
ticks
30.0

BUTTON
10
10
73
43
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
73
10
136
43
GO!
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
135
10
198
43
Step
go\n;ask turtles [\n;  set speak? TRUE\n;  getColor\n;]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
10
42
102
75
N-agents
N-agents
10
1000
129.0
1
1
NIL
HORIZONTAL

SLIDER
10
352
102
385
n-neis
n-neis
1
500
64.0
1
1
NIL
HORIZONTAL

SLIDER
102
352
194
385
p-random
p-random
0
0.5
0.27
0.01
1
NIL
HORIZONTAL

SLIDER
102
42
194
75
opinions
opinions
1
50
4.0
1
1
NIL
HORIZONTAL

CHOOSER
10
108
102
153
model
model
"HK"
0

SLIDER
10
274
138
307
p-speaking-level
p-speaking-level
0
1
0.543
0.001
1
NIL
HORIZONTAL

SLIDER
10
196
102
229
boundary
boundary
0.01
1
0.22
0.01
1
NIL
HORIZONTAL

INPUTBOX
674
10
781
70
RS
24.0
1
0
Number

SWITCH
781
10
874
43
set-seed?
set-seed?
0
1
-1000

BUTTON
435
491
490
524
inspect
inspect turtle 0\nask turtle 0 [print opinion-position]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
490
491
545
524
sizes
show sort remove-duplicates [count turtles-here] of turtles \nask patches [set plabel count turtles-here]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
379
491
434
524
HIDE!
\nask turtles [set hidden? TRUE]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
379
458
434
491
SHOW!
\nask turtles [set hidden? FALSE]\nask turtles [set size (max-pxcor / 10)]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

CHOOSER
10
229
102
274
boundary-drawn
boundary-drawn
"constant" "uniform"
1

PLOT
967
256
1127
376
Distribution of 'Uncertainty'
NIL
NIL
0.0
1.0
0.0
10.0
true
false
"" ""
PENS
"default" 0.05 1 -16777216 true "" "histogram [Uncertainty] of agents"

BUTTON
434
458
522
491
avg. Unretainty
show mean [Uncertainty] of turtles
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
521
458
583
491
Hide links
ask links [set hidden? TRUE]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
583
458
643
491
Show links
ask comms [set hidden? FALSE]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
931
10
1023
43
X-opinion
X-opinion
1
50
1.0
1
1
NIL
HORIZONTAL

SLIDER
931
42
1023
75
Y-opinion
Y-opinion
1
50
4.0
1
1
NIL
HORIZONTAL

PLOT
674
224
968
432
Developement of opinions
NIL
NIL
0.0
10.0
-0.1
0.1
true
true
"" ""
PENS
"Op01" 1.0 0 -16777216 true "" "plot mean [item 0 opinion-position] of agents"
"Op02" 1.0 0 -7500403 true "" "if opinions >= 2 [plot mean [item 1 opinion-position] of agents]"
"Op03" 1.0 0 -2674135 true "" "if opinions >= 3 [plot mean [item 2 opinion-position] of agents]"
"Op04" 1.0 0 -955883 true "" "if opinions >= 4 [plot mean [item 3 opinion-position] of agents]"
"Op05" 1.0 0 -6459832 true "" "if opinions >= 5 [plot mean [item 4 opinion-position] of turtles]"
"Op06" 1.0 0 -1184463 true "" "if opinions >= 6 [plot mean [item 5 opinion-position] of turtles]"
"Op07" 1.0 0 -10899396 true "" "if opinions >= 7 [plot mean [item 6 opinion-position] of turtles]"
"Op08" 1.0 0 -13840069 true "" "if opinions >= 8 [plot mean [item 7 opinion-position] of turtles]"
"Op09" 1.0 0 -14835848 true "" "if opinions >= 9 [plot mean [item 8 opinion-position] of turtles]"
"Op10" 1.0 0 -11221820 true "" "if opinions >= 10 [plot mean [item 9 opinion-position] of turtles]"
"Op11" 1.0 0 -13791810 true "" "if opinions >= 11 [plot mean [item 10 opinion-position] of turtles]"

SLIDER
931
74
1023
107
updating
updating
1
50
1.0
1
1
NIL
HORIZONTAL

PLOT
674
107
998
227
Stability of turtles (average)
NIL
NIL
0.0
10.0
0.0
0.3
true
true
"" ""
PENS
"Turtles" 1.0 0 -16777216 true "" "plot mean [mean Record] of agents"
"Main record" 1.0 0 -2674135 true "" "plot mean main-Record"
"Normalized" 1.0 0 -13791810 true "" "plot normalized_polarisation"
"ESBG" 1.0 0 -11221820 true "" "plot ESBG_polarisation"

BUTTON
1403
161
1475
194
avg. Record
show mean [mean Record] of agents
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

MONITOR
1404
117
1476
162
avg. Record
mean [mean Record] of agents
3
1
11

SLIDER
1023
10
1134
43
record-length
record-length
10
100
40.0
1
1
NIL
HORIZONTAL

CHOOSER
10
152
102
197
mode
mode
"openly-listen" "vaguely-speak"
0

PLOT
967
375
1127
495
Distribution of 'Tolerance'
NIL
NIL
0.0
1.0
0.0
10.0
true
false
"" ""
PENS
"default" 0.05 1 -16777216 true "" "histogram [Tolerance] of agents"

INPUTBOX
1171
10
1423
116
file-name-core
1_129_0.1_16_2_2_0.4_uniform_1_uniform_openly-listen
1
0
String

BUTTON
1319
116
1404
149
RECORD!
record-state-of-simulation
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SWITCH
1171
116
1319
149
construct-name?
construct-name?
1
1
-1000

SWITCH
1171
148
1274
181
record?
record?
1
1
-1000

SLIDER
1023
42
1134
75
max-ticks
max-ticks
100
10000
1000.0
100
1
NIL
HORIZONTAL

SWITCH
1274
148
1396
181
HK-benchmark?
HK-benchmark?
1
1
-1000

CHOOSER
10
307
107
352
p-speaking-drawn
p-speaking-drawn
"constant" "uniform" "function"
1

PLOT
1127
375
1287
495
Distribution of 'Outspokeness'
NIL
NIL
0.0
1.0
0.0
10.0
true
false
"" ""
PENS
"default" 0.05 1 -16777216 true "" "histogram [P-speaking] of agents"

SLIDER
1023
74
1171
107
record-each-n-steps
record-each-n-steps
100
10000
1000.0
100
1
NIL
HORIZONTAL

SWITCH
781
43
929
76
avoid-redundancies?
avoid-redundancies?
1
1
-1000

SLIDER
10
75
126
108
tolerance-level
tolerance-level
0
1.1
1.08
0.01
1
NIL
HORIZONTAL

CHOOSER
126
75
218
120
tolerance-drawn
tolerance-drawn
"constant" "uniform"
1

PLOT
767
432
967
552
Number of network changes
NIL
NIL
0.0
10.0
0.0
100.0
true
true
"" ""
PENS
"changed" 1.0 0 -16777216 true "" "plot network-changes / count agents * 100"
"satisfied" 1.0 0 -2674135 true "" "plot count (agents with [Satisfied?]) / count agents * 100"

SLIDER
10
429
124
462
conformity-level
conformity-level
0
1
0.36
0.01
1
NIL
HORIZONTAL

CHOOSER
10
385
129
430
conformity-drawn
conformity-drawn
"constant" "uniform"
1

CHOOSER
128
385
220
430
network-change
network-change
"link" "community"
0

MONITOR
674
432
768
477
NIL
network-changes
17
1
11

SWITCH
10
462
150
495
cut-links-randomly?
cut-links-randomly?
1
1
-1000

BUTTON
545
491
620
524
polarisation
compute-polarisation-repeatedly
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

INPUTBOX
1316
182
1389
242
N_centroids
2.0
1
0
Number

SLIDER
1021
193
1171
226
Centroids_change
Centroids_change
0.00000001
0.001
1.0E-5
0.00001
1
NIL
HORIZONTAL

PLOT
1127
256
1287
376
Distribution of 'Conformity'
NIL
NIL
0.0
1.0
0.0
10.0
true
false
"" ""
PENS
"default" 0.05 1 -16777216 true "" "histogram [Conformity] of agents"

MONITOR
643
507
768
552
NIL
normalized_polarisation
17
1
11

MONITOR
967
507
1068
552
NIL
ESBG_polarisation
17
1
11

SWITCH
1021
225
1171
258
centroid_color?
centroid_color?
1
1
-1000

SWITCH
1171
224
1319
257
killing_centroids?
killing_centroids?
0
1
-1000

SLIDER
102
196
206
229
id_threshold
id_threshold
0.01
1
0.49
0.01
1
NIL
HORIZONTAL

SLIDER
1001
127
1171
160
polarisation-each-n-steps
polarisation-each-n-steps
0
10000
1010.0
10
1
NIL
HORIZONTAL

SLIDER
1171
191
1276
224
polar_repeats
polar_repeats
1
100
20.0
1
1
NIL
HORIZONTAL

SWITCH
102
163
217
196
use_identity?
use_identity?
0
1
-1000

SLIDER
1021
160
1172
193
d_threshold
d_threshold
0
1
0.8
.01
1
NIL
HORIZONTAL

SLIDER
643
476
767
509
ESBG_furthest_out
ESBG_furthest_out
0
100
5.0
1
1
NIL
HORIZONTAL

CHOOSER
102
229
194
274
threshold_drawn
threshold_drawn
"constant" "uniform"
1

PLOT
1287
256
1447
378
Distribution of 'Group threshold'
NIL
NIL
0.0
1.0
0.0
10.0
true
false
"" ""
PENS
"default" 0.05 1 -16777216 true "" "histogram [group-threshold] of agents"

PLOT
1287
375
1447
495
Distribution od 'Size of identity group'
NIL
NIL
0.0
130.0
0.0
10.0
true
false
"" ""
PENS
"default" 10.0 1 -16777216 true "" "histogram [count Identity-group] of agents"

CHOOSER
107
307
199
352
identity_type
identity_type
"global" "individual"
0

SWITCH
149
462
306
495
create-links-randomly?
create-links-randomly?
1
1
-1000

SLIDER
118
429
219
462
min-comm-neis
min-comm-neis
0
10
5.0
1
1
NIL
HORIZONTAL

SLIDER
10
494
192
527
dissatisfied_updates_opinion
dissatisfied_updates_opinion
0
1
0.41
0.01
1
NIL
HORIZONTAL

SWITCH
191
494
345
527
use_opponents_ratio?
use_opponents_ratio?
1
1
-1000

SLIDER
102
130
217
163
identity_levels
identity_levels
2
10
2.0
1
1
NIL
HORIZONTAL

CHOOSER
191
526
329
571
draw_id_threshold
draw_id_threshold
"uniform" "left-skewed" "right-skewed"
0

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)  

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.2.2
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="secondTry" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>ticks</metric>
    <enumeratedValueSet variable="max-ticks">
      <value value="5000"/>
    </enumeratedValueSet>
    <steppedValueSet variable="RS" first="1" step="1" last="10"/>
    <enumeratedValueSet variable="N-agents">
      <value value="129"/>
      <value value="257"/>
      <value value="513"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-random">
      <value value="0.05"/>
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-neis">
      <value value="4"/>
      <value value="16"/>
      <value value="64"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="opinions">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="updating">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary">
      <value value="0.1"/>
      <value value="0.2"/>
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary-drawn">
      <value value="&quot;constant&quot;"/>
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking">
      <value value="0.5"/>
      <value value="0.9"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mode">
      <value value="&quot;openly-listen&quot;"/>
      <value value="&quot;vaguely-speak&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="smallest-component">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X-opinion">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Y-opinion">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="&quot;HK&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-boundary">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sigma">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="construct-name?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mu">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-length">
      <value value="100"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="thirdTry" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>ticks</metric>
    <enumeratedValueSet variable="max-ticks">
      <value value="5000"/>
    </enumeratedValueSet>
    <steppedValueSet variable="RS" first="1" step="1" last="10"/>
    <enumeratedValueSet variable="N-agents">
      <value value="129"/>
      <value value="257"/>
      <value value="513"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-random">
      <value value="0.05"/>
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-neis">
      <value value="4"/>
      <value value="16"/>
      <value value="64"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="opinions">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="updating">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary">
      <value value="0.1"/>
      <value value="0.2"/>
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary-drawn">
      <value value="&quot;constant&quot;"/>
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking">
      <value value="0.5"/>
      <value value="0.9"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mode">
      <value value="&quot;openly-listen&quot;"/>
      <value value="&quot;vaguely-speak&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="smallest-component">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X-opinion">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Y-opinion">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="&quot;HK&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-boundary">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sigma">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="construct-name?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mu">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-length">
      <value value="100"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="HK-benchmark" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>ticks</metric>
    <enumeratedValueSet variable="record?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="construct-name?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="5000"/>
    </enumeratedValueSet>
    <steppedValueSet variable="RS" first="1" step="1" last="10"/>
    <enumeratedValueSet variable="N-agents">
      <value value="129"/>
      <value value="257"/>
      <value value="513"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-random">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="HK-benchmark?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="opinions">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="updating">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary">
      <value value="0.1"/>
      <value value="0.2"/>
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary-drawn">
      <value value="&quot;constant&quot;"/>
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mode">
      <value value="&quot;openly-listen&quot;"/>
      <value value="&quot;vaguely-speak&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X-opinion">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Y-opinion">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="&quot;HK&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-boundary">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sigma">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mu">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-length">
      <value value="100"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="HK-benchmark-SoS" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>ticks</metric>
    <enumeratedValueSet variable="record?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="construct-name?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="5000"/>
    </enumeratedValueSet>
    <steppedValueSet variable="RS" first="1" step="1" last="10"/>
    <enumeratedValueSet variable="N-agents">
      <value value="129"/>
      <value value="257"/>
      <value value="513"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-random">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="HK-benchmark?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="opinions">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="updating">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary">
      <value value="0.1"/>
      <value value="0.2"/>
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary-drawn">
      <value value="&quot;constant&quot;"/>
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking">
      <value value="0.9"/>
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mode">
      <value value="&quot;openly-listen&quot;"/>
      <value value="&quot;vaguely-speak&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X-opinion">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Y-opinion">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="&quot;HK&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-boundary">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sigma">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mu">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-length">
      <value value="100"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="HK-benchmarkV02" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>ticks</metric>
    <enumeratedValueSet variable="record?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="construct-name?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="5000"/>
    </enumeratedValueSet>
    <steppedValueSet variable="RS" first="11" step="1" last="73"/>
    <enumeratedValueSet variable="N-agents">
      <value value="129"/>
      <value value="257"/>
      <value value="513"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-random">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="HK-benchmark?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="opinions">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="updating">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary">
      <value value="0.1"/>
      <value value="0.2"/>
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary-drawn">
      <value value="&quot;constant&quot;"/>
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking">
      <value value="1"/>
      <value value="0.9"/>
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mode">
      <value value="&quot;openly-listen&quot;"/>
      <value value="&quot;vaguely-speak&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X-opinion">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Y-opinion">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="&quot;HK&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-boundary">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sigma">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mu">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-length">
      <value value="100"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="NewExperiment" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="1001"/>
    <enumeratedValueSet variable="record?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="HK-benchmark?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="&quot;HK&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-each-n-steps">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="construct-name?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-length">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="avoid-redundancies?">
      <value value="true"/>
    </enumeratedValueSet>
    <steppedValueSet variable="RS" first="101" step="1" last="110"/>
    <enumeratedValueSet variable="N-agents">
      <value value="101"/>
      <value value="1001"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-random">
      <value value="0"/>
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-neis">
      <value value="5"/>
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="opinions">
      <value value="1"/>
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="updating">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary">
      <value value="0.1"/>
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary-drawn">
      <value value="&quot;constant&quot;"/>
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-boundary">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking-level">
      <value value="0.5"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking-drawn">
      <value value="&quot;constant&quot;"/>
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mode">
      <value value="&quot;openly-listen&quot;"/>
      <value value="&quot;vaguely-speak&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X-opinion">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Y-opinion">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mu">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sigma">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bias-margin">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bias-of-bias">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bias-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="social-bias">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="regularExperiment" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="366"/>
    <metric>betweenness_start</metric>
    <metric>eigenvector_start</metric>
    <metric>clustering_start</metric>
    <metric>mean_path_start</metric>
    <metric>normalized_polarization_start</metric>
    <metric>ESBSG_polarization_start</metric>
    <metric>mean_op1_start</metric>
    <metric>mean_op2_start</metric>
    <metric>sd_op1_start</metric>
    <metric>sd_op2_start</metric>
    <metric>median_op1_start</metric>
    <metric>median_op2_start</metric>
    <metric>lower_op1_start</metric>
    <metric>lower_op2_start</metric>
    <metric>upper_op1_start</metric>
    <metric>upper_op2_start</metric>
    <metric>betweenness_final</metric>
    <metric>eigenvector_final</metric>
    <metric>clustering_final</metric>
    <metric>mean_path_final</metric>
    <metric>normalized_polarization_final</metric>
    <metric>ESBSG_polarization_final</metric>
    <metric>mean_op1_final</metric>
    <metric>mean_op2_final</metric>
    <metric>sd_op1_final</metric>
    <metric>sd_op2_final</metric>
    <metric>median_op1_final</metric>
    <metric>median_op2_final</metric>
    <metric>lower_op1_final</metric>
    <metric>lower_op2_final</metric>
    <metric>upper_op1_final</metric>
    <metric>upper_op2_final</metric>
    <steppedValueSet variable="RS" first="1" step="1" last="30"/>
    <enumeratedValueSet variable="cut-links-randomly?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="opinions">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="updating">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary">
      <value value="0.22"/>
      <value value="0.28"/>
      <value value="0.34"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="polar_repeats">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Y-opinion">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="draw_id_threshold">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="&quot;HK&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="365"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Centroids_change">
      <value value="1.0E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-random">
      <value value="0.27"/>
      <value value="0.34"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tolerance-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-neis">
      <value value="64"/>
      <value value="52"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="use_opponents_ratio?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="polarisation-each-n-steps">
      <value value="1010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mode">
      <value value="&quot;openly-listen&quot;"/>
      <value value="&quot;vaguely-speak&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="d_threshold">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-each-n-steps">
      <value value="365"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ESBG_furthest_out">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="identity_levels">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-length">
      <value value="40"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="id_threshold">
      <value value="0.39"/>
      <value value="0.49"/>
      <value value="0.59"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="network-change">
      <value value="&quot;link&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tolerance-level">
      <value value="1.08"/>
      <value value="0.648"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking-level">
      <value value="0.43"/>
      <value value="0.543"/>
      <value value="0.65"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="centroid_color?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="use_identity?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conformity-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N_centroids">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dissatisfied_updates_opinion">
      <value value="0.41"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="create-links-randomly?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="killing_centroids?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N-agents">
      <value value="129"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="construct-name?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="avoid-redundancies?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="HK-benchmark?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-comm-neis">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conformity-level">
      <value value="0.36"/>
      <value value="0.45"/>
      <value value="0.54"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold_drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="identity_type">
      <value value="&quot;global&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X-opinion">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="regularExperimentPart02" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="366"/>
    <metric>betweenness_start</metric>
    <metric>eigenvector_start</metric>
    <metric>clustering_start</metric>
    <metric>mean_path_start</metric>
    <metric>normalized_polarization_start</metric>
    <metric>ESBSG_polarization_start</metric>
    <metric>mean_op1_start</metric>
    <metric>mean_op2_start</metric>
    <metric>sd_op1_start</metric>
    <metric>sd_op2_start</metric>
    <metric>median_op1_start</metric>
    <metric>median_op2_start</metric>
    <metric>lower_op1_start</metric>
    <metric>lower_op2_start</metric>
    <metric>upper_op1_start</metric>
    <metric>upper_op2_start</metric>
    <metric>betweenness_final</metric>
    <metric>eigenvector_final</metric>
    <metric>clustering_final</metric>
    <metric>mean_path_final</metric>
    <metric>normalized_polarization_final</metric>
    <metric>ESBSG_polarization_final</metric>
    <metric>mean_op1_final</metric>
    <metric>mean_op2_final</metric>
    <metric>sd_op1_final</metric>
    <metric>sd_op2_final</metric>
    <metric>median_op1_final</metric>
    <metric>median_op2_final</metric>
    <metric>lower_op1_final</metric>
    <metric>lower_op2_final</metric>
    <metric>upper_op1_final</metric>
    <metric>upper_op2_final</metric>
    <steppedValueSet variable="RS" first="31" step="1" last="60"/>
    <enumeratedValueSet variable="cut-links-randomly?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="opinions">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="updating">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary">
      <value value="0.22"/>
      <value value="0.28"/>
      <value value="0.34"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="polar_repeats">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Y-opinion">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="draw_id_threshold">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="&quot;HK&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="365"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Centroids_change">
      <value value="1.0E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-random">
      <value value="0.27"/>
      <value value="0.34"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tolerance-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-neis">
      <value value="64"/>
      <value value="52"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="use_opponents_ratio?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="polarisation-each-n-steps">
      <value value="1010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mode">
      <value value="&quot;openly-listen&quot;"/>
      <value value="&quot;vaguely-speak&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="d_threshold">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-each-n-steps">
      <value value="365"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ESBG_furthest_out">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="identity_levels">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-length">
      <value value="40"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="id_threshold">
      <value value="0.39"/>
      <value value="0.49"/>
      <value value="0.59"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="network-change">
      <value value="&quot;link&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tolerance-level">
      <value value="1.08"/>
      <value value="0.648"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking-level">
      <value value="0.43"/>
      <value value="0.543"/>
      <value value="0.65"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="centroid_color?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="use_identity?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conformity-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N_centroids">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dissatisfied_updates_opinion">
      <value value="0.41"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="create-links-randomly?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="killing_centroids?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N-agents">
      <value value="129"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="construct-name?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="avoid-redundancies?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="HK-benchmark?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-comm-neis">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conformity-level">
      <value value="0.36"/>
      <value value="0.45"/>
      <value value="0.54"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold_drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="identity_type">
      <value value="&quot;global&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X-opinion">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="regularExperimentPart03" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="366"/>
    <metric>betweenness_start</metric>
    <metric>eigenvector_start</metric>
    <metric>clustering_start</metric>
    <metric>mean_path_start</metric>
    <metric>normalized_polarization_start</metric>
    <metric>ESBSG_polarization_start</metric>
    <metric>mean_op1_start</metric>
    <metric>mean_op2_start</metric>
    <metric>sd_op1_start</metric>
    <metric>sd_op2_start</metric>
    <metric>median_op1_start</metric>
    <metric>median_op2_start</metric>
    <metric>lower_op1_start</metric>
    <metric>lower_op2_start</metric>
    <metric>upper_op1_start</metric>
    <metric>upper_op2_start</metric>
    <metric>betweenness_final</metric>
    <metric>eigenvector_final</metric>
    <metric>clustering_final</metric>
    <metric>mean_path_final</metric>
    <metric>normalized_polarization_final</metric>
    <metric>ESBSG_polarization_final</metric>
    <metric>mean_op1_final</metric>
    <metric>mean_op2_final</metric>
    <metric>sd_op1_final</metric>
    <metric>sd_op2_final</metric>
    <metric>median_op1_final</metric>
    <metric>median_op2_final</metric>
    <metric>lower_op1_final</metric>
    <metric>lower_op2_final</metric>
    <metric>upper_op1_final</metric>
    <metric>upper_op2_final</metric>
    <steppedValueSet variable="RS" first="1" step="1" last="60"/>
    <enumeratedValueSet variable="cut-links-randomly?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="opinions">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="updating">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary">
      <value value="0.22"/>
      <value value="0.28"/>
      <value value="0.34"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="polar_repeats">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Y-opinion">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="draw_id_threshold">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="&quot;HK&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="365"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Centroids_change">
      <value value="1.0E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-random">
      <value value="0.27"/>
      <value value="0.34"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tolerance-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-neis">
      <value value="64"/>
      <value value="52"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="use_opponents_ratio?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="polarisation-each-n-steps">
      <value value="1010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mode">
      <value value="&quot;openly-listen&quot;"/>
      <value value="&quot;vaguely-speak&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="d_threshold">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-each-n-steps">
      <value value="365"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ESBG_furthest_out">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="identity_levels">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-length">
      <value value="40"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="id_threshold">
      <value value="0.49"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="network-change">
      <value value="&quot;link&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tolerance-level">
      <value value="1.08"/>
      <value value="0.648"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking-level">
      <value value="0.43"/>
      <value value="0.543"/>
      <value value="0.65"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="centroid_color?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="use_identity?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conformity-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N_centroids">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dissatisfied_updates_opinion">
      <value value="0.41"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="create-links-randomly?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="killing_centroids?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N-agents">
      <value value="129"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="construct-name?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="avoid-redundancies?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="HK-benchmark?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-comm-neis">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conformity-level">
      <value value="0.36"/>
      <value value="0.45"/>
      <value value="0.54"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold_drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="identity_type">
      <value value="&quot;global&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X-opinion">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="regularExperimentPart03LONG" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="3651"/>
    <metric>betweenness_start</metric>
    <metric>eigenvector_start</metric>
    <metric>clustering_start</metric>
    <metric>mean_path_start</metric>
    <metric>normalized_polarization_start</metric>
    <metric>ESBSG_polarization_start</metric>
    <metric>mean_op1_start</metric>
    <metric>mean_op2_start</metric>
    <metric>sd_op1_start</metric>
    <metric>sd_op2_start</metric>
    <metric>median_op1_start</metric>
    <metric>median_op2_start</metric>
    <metric>lower_op1_start</metric>
    <metric>lower_op2_start</metric>
    <metric>upper_op1_start</metric>
    <metric>upper_op2_start</metric>
    <metric>betweenness_final</metric>
    <metric>eigenvector_final</metric>
    <metric>clustering_final</metric>
    <metric>mean_path_final</metric>
    <metric>normalized_polarization_final</metric>
    <metric>ESBSG_polarization_final</metric>
    <metric>mean_op1_final</metric>
    <metric>mean_op2_final</metric>
    <metric>sd_op1_final</metric>
    <metric>sd_op2_final</metric>
    <metric>median_op1_final</metric>
    <metric>median_op2_final</metric>
    <metric>lower_op1_final</metric>
    <metric>lower_op2_final</metric>
    <metric>upper_op1_final</metric>
    <metric>upper_op2_final</metric>
    <steppedValueSet variable="RS" first="1" step="1" last="60"/>
    <enumeratedValueSet variable="cut-links-randomly?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="opinions">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="updating">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary">
      <value value="0.22"/>
      <value value="0.28"/>
      <value value="0.34"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="polar_repeats">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Y-opinion">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="draw_id_threshold">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="&quot;HK&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="3650"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Centroids_change">
      <value value="1.0E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-random">
      <value value="0.27"/>
      <value value="0.34"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tolerance-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-neis">
      <value value="64"/>
      <value value="52"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="use_opponents_ratio?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="polarisation-each-n-steps">
      <value value="10010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mode">
      <value value="&quot;openly-listen&quot;"/>
      <value value="&quot;vaguely-speak&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="d_threshold">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-each-n-steps">
      <value value="3650"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ESBG_furthest_out">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="identity_levels">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-length">
      <value value="40"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="id_threshold">
      <value value="0.49"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="network-change">
      <value value="&quot;link&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tolerance-level">
      <value value="1.08"/>
      <value value="0.648"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking-level">
      <value value="0.43"/>
      <value value="0.543"/>
      <value value="0.65"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="centroid_color?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="use_identity?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conformity-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N_centroids">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dissatisfied_updates_opinion">
      <value value="0.41"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="create-links-randomly?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="killing_centroids?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N-agents">
      <value value="129"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="construct-name?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="avoid-redundancies?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="HK-benchmark?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-comm-neis">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conformity-level">
      <value value="0.36"/>
      <value value="0.45"/>
      <value value="0.54"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold_drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="identity_type">
      <value value="&quot;global&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X-opinion">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="regularExperimentPart04" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="366"/>
    <metric>betweenness_start</metric>
    <metric>eigenvector_start</metric>
    <metric>clustering_start</metric>
    <metric>mean_path_start</metric>
    <metric>normalized_polarization_start</metric>
    <metric>ESBSG_polarization_start</metric>
    <metric>mean_op1_start</metric>
    <metric>mean_op2_start</metric>
    <metric>sd_op1_start</metric>
    <metric>sd_op2_start</metric>
    <metric>median_op1_start</metric>
    <metric>median_op2_start</metric>
    <metric>lower_op1_start</metric>
    <metric>lower_op2_start</metric>
    <metric>upper_op1_start</metric>
    <metric>upper_op2_start</metric>
    <metric>betweenness_final</metric>
    <metric>eigenvector_final</metric>
    <metric>clustering_final</metric>
    <metric>mean_path_final</metric>
    <metric>normalized_polarization_final</metric>
    <metric>ESBSG_polarization_final</metric>
    <metric>mean_op1_final</metric>
    <metric>mean_op2_final</metric>
    <metric>sd_op1_final</metric>
    <metric>sd_op2_final</metric>
    <metric>median_op1_final</metric>
    <metric>median_op2_final</metric>
    <metric>lower_op1_final</metric>
    <metric>lower_op2_final</metric>
    <metric>upper_op1_final</metric>
    <metric>upper_op2_final</metric>
    <steppedValueSet variable="RS" first="1" step="1" last="60"/>
    <enumeratedValueSet variable="n-neis">
      <value value="52"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tolerance-level">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cut-links-randomly?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="opinions">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="updating">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary">
      <value value="0.22"/>
      <value value="0.28"/>
      <value value="0.34"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="polar_repeats">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Y-opinion">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="draw_id_threshold">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="&quot;HK&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="365"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Centroids_change">
      <value value="1.0E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-random">
      <value value="0.27"/>
      <value value="0.34"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tolerance-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="use_opponents_ratio?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="polarisation-each-n-steps">
      <value value="1010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mode">
      <value value="&quot;openly-listen&quot;"/>
      <value value="&quot;vaguely-speak&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="d_threshold">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-each-n-steps">
      <value value="365"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ESBG_furthest_out">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="identity_levels">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-length">
      <value value="40"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="id_threshold">
      <value value="0.39"/>
      <value value="0.49"/>
      <value value="0.59"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="network-change">
      <value value="&quot;link&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking-level">
      <value value="0.43"/>
      <value value="0.543"/>
      <value value="0.65"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="centroid_color?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="use_identity?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conformity-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N_centroids">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dissatisfied_updates_opinion">
      <value value="0.41"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="create-links-randomly?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="killing_centroids?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N-agents">
      <value value="129"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="construct-name?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="avoid-redundancies?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="HK-benchmark?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-comm-neis">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conformity-level">
      <value value="0.36"/>
      <value value="0.45"/>
      <value value="0.54"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold_drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="identity_type">
      <value value="&quot;global&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X-opinion">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="regularExperimentPart06" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="366"/>
    <metric>betweenness_start</metric>
    <metric>eigenvector_start</metric>
    <metric>clustering_start</metric>
    <metric>mean_path_start</metric>
    <metric>normalized_polarization_start</metric>
    <metric>ESBSG_polarization_start</metric>
    <metric>mean_op1_start</metric>
    <metric>mean_op2_start</metric>
    <metric>sd_op1_start</metric>
    <metric>sd_op2_start</metric>
    <metric>median_op1_start</metric>
    <metric>median_op2_start</metric>
    <metric>lower_op1_start</metric>
    <metric>lower_op2_start</metric>
    <metric>upper_op1_start</metric>
    <metric>upper_op2_start</metric>
    <metric>betweenness_final</metric>
    <metric>eigenvector_final</metric>
    <metric>clustering_final</metric>
    <metric>mean_path_final</metric>
    <metric>normalized_polarization_final</metric>
    <metric>ESBSG_polarization_final</metric>
    <metric>mean_op1_final</metric>
    <metric>mean_op2_final</metric>
    <metric>sd_op1_final</metric>
    <metric>sd_op2_final</metric>
    <metric>median_op1_final</metric>
    <metric>median_op2_final</metric>
    <metric>lower_op1_final</metric>
    <metric>lower_op2_final</metric>
    <metric>upper_op1_final</metric>
    <metric>upper_op2_final</metric>
    <steppedValueSet variable="RS" first="1" step="1" last="60"/>
    <enumeratedValueSet variable="n-neis">
      <value value="13"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tolerance-level">
      <value value="0.648"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cut-links-randomly?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="opinions">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="updating">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary">
      <value value="0.22"/>
      <value value="0.28"/>
      <value value="0.34"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="polar_repeats">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Y-opinion">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="draw_id_threshold">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="&quot;HK&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="365"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Centroids_change">
      <value value="1.0E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-random">
      <value value="0.27"/>
      <value value="0.34"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tolerance-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="use_opponents_ratio?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="polarisation-each-n-steps">
      <value value="1010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mode">
      <value value="&quot;openly-listen&quot;"/>
      <value value="&quot;vaguely-speak&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="d_threshold">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-each-n-steps">
      <value value="365"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ESBG_furthest_out">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="identity_levels">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-length">
      <value value="40"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="id_threshold">
      <value value="0.39"/>
      <value value="0.49"/>
      <value value="0.59"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="network-change">
      <value value="&quot;link&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking-level">
      <value value="0.43"/>
      <value value="0.543"/>
      <value value="0.65"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="centroid_color?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="use_identity?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conformity-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N_centroids">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dissatisfied_updates_opinion">
      <value value="0.41"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="create-links-randomly?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="killing_centroids?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N-agents">
      <value value="129"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="construct-name?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="avoid-redundancies?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="HK-benchmark?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-comm-neis">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conformity-level">
      <value value="0.36"/>
      <value value="0.45"/>
      <value value="0.54"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold_drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="identity_type">
      <value value="&quot;global&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X-opinion">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="regularExperimentPart05" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="366"/>
    <metric>betweenness_start</metric>
    <metric>eigenvector_start</metric>
    <metric>clustering_start</metric>
    <metric>mean_path_start</metric>
    <metric>normalized_polarization_start</metric>
    <metric>ESBSG_polarization_start</metric>
    <metric>mean_op1_start</metric>
    <metric>mean_op2_start</metric>
    <metric>sd_op1_start</metric>
    <metric>sd_op2_start</metric>
    <metric>median_op1_start</metric>
    <metric>median_op2_start</metric>
    <metric>lower_op1_start</metric>
    <metric>lower_op2_start</metric>
    <metric>upper_op1_start</metric>
    <metric>upper_op2_start</metric>
    <metric>betweenness_final</metric>
    <metric>eigenvector_final</metric>
    <metric>clustering_final</metric>
    <metric>mean_path_final</metric>
    <metric>normalized_polarization_final</metric>
    <metric>ESBSG_polarization_final</metric>
    <metric>mean_op1_final</metric>
    <metric>mean_op2_final</metric>
    <metric>sd_op1_final</metric>
    <metric>sd_op2_final</metric>
    <metric>median_op1_final</metric>
    <metric>median_op2_final</metric>
    <metric>lower_op1_final</metric>
    <metric>lower_op2_final</metric>
    <metric>upper_op1_final</metric>
    <metric>upper_op2_final</metric>
    <steppedValueSet variable="RS" first="1" step="1" last="60"/>
    <enumeratedValueSet variable="n-neis">
      <value value="13"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tolerance-level">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cut-links-randomly?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="opinions">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="updating">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary">
      <value value="0.22"/>
      <value value="0.28"/>
      <value value="0.34"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="polar_repeats">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Y-opinion">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="draw_id_threshold">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="&quot;HK&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="365"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Centroids_change">
      <value value="1.0E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-random">
      <value value="0.27"/>
      <value value="0.34"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tolerance-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="use_opponents_ratio?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="polarisation-each-n-steps">
      <value value="1010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mode">
      <value value="&quot;openly-listen&quot;"/>
      <value value="&quot;vaguely-speak&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="d_threshold">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-each-n-steps">
      <value value="365"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ESBG_furthest_out">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="identity_levels">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-length">
      <value value="40"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="id_threshold">
      <value value="0.39"/>
      <value value="0.49"/>
      <value value="0.59"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="network-change">
      <value value="&quot;link&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking-level">
      <value value="0.43"/>
      <value value="0.543"/>
      <value value="0.65"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="centroid_color?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="use_identity?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conformity-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N_centroids">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dissatisfied_updates_opinion">
      <value value="0.41"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="create-links-randomly?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="killing_centroids?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N-agents">
      <value value="129"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="construct-name?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="avoid-redundancies?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="HK-benchmark?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-comm-neis">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conformity-level">
      <value value="0.36"/>
      <value value="0.45"/>
      <value value="0.54"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold_drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="identity_type">
      <value value="&quot;global&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X-opinion">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="regularExperimentPart01LONG01" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="3651"/>
    <metric>betweenness_start</metric>
    <metric>eigenvector_start</metric>
    <metric>clustering_start</metric>
    <metric>mean_path_start</metric>
    <metric>normalized_polarization_start</metric>
    <metric>ESBSG_polarization_start</metric>
    <metric>mean_op1_start</metric>
    <metric>mean_op2_start</metric>
    <metric>sd_op1_start</metric>
    <metric>sd_op2_start</metric>
    <metric>median_op1_start</metric>
    <metric>median_op2_start</metric>
    <metric>lower_op1_start</metric>
    <metric>lower_op2_start</metric>
    <metric>upper_op1_start</metric>
    <metric>upper_op2_start</metric>
    <metric>betweenness_final</metric>
    <metric>eigenvector_final</metric>
    <metric>clustering_final</metric>
    <metric>mean_path_final</metric>
    <metric>normalized_polarization_final</metric>
    <metric>ESBSG_polarization_final</metric>
    <metric>mean_op1_final</metric>
    <metric>mean_op2_final</metric>
    <metric>sd_op1_final</metric>
    <metric>sd_op2_final</metric>
    <metric>median_op1_final</metric>
    <metric>median_op2_final</metric>
    <metric>lower_op1_final</metric>
    <metric>lower_op2_final</metric>
    <metric>upper_op1_final</metric>
    <metric>upper_op2_final</metric>
    <steppedValueSet variable="RS" first="1" step="1" last="3"/>
    <enumeratedValueSet variable="cut-links-randomly?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="opinions">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="updating">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary">
      <value value="0.22"/>
      <value value="0.28"/>
      <value value="0.34"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="polar_repeats">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Y-opinion">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="draw_id_threshold">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="&quot;HK&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="3650"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Centroids_change">
      <value value="1.0E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-random">
      <value value="0.27"/>
      <value value="0.34"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tolerance-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-neis">
      <value value="64"/>
      <value value="52"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="use_opponents_ratio?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="polarisation-each-n-steps">
      <value value="10010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mode">
      <value value="&quot;openly-listen&quot;"/>
      <value value="&quot;vaguely-speak&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="d_threshold">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-each-n-steps">
      <value value="3650"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ESBG_furthest_out">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="identity_levels">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-length">
      <value value="40"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="id_threshold">
      <value value="0.39"/>
      <value value="0.49"/>
      <value value="0.59"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="network-change">
      <value value="&quot;link&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tolerance-level">
      <value value="1.08"/>
      <value value="0.648"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking-level">
      <value value="0.43"/>
      <value value="0.543"/>
      <value value="0.65"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="centroid_color?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="use_identity?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conformity-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N_centroids">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dissatisfied_updates_opinion">
      <value value="0.41"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="create-links-randomly?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="killing_centroids?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N-agents">
      <value value="129"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="construct-name?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="avoid-redundancies?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="HK-benchmark?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-comm-neis">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conformity-level">
      <value value="0.36"/>
      <value value="0.45"/>
      <value value="0.54"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold_drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="identity_type">
      <value value="&quot;global&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X-opinion">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="regularExperimentPart01LONG02" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="3651"/>
    <metric>betweenness_start</metric>
    <metric>eigenvector_start</metric>
    <metric>clustering_start</metric>
    <metric>mean_path_start</metric>
    <metric>normalized_polarization_start</metric>
    <metric>ESBSG_polarization_start</metric>
    <metric>mean_op1_start</metric>
    <metric>mean_op2_start</metric>
    <metric>sd_op1_start</metric>
    <metric>sd_op2_start</metric>
    <metric>median_op1_start</metric>
    <metric>median_op2_start</metric>
    <metric>lower_op1_start</metric>
    <metric>lower_op2_start</metric>
    <metric>upper_op1_start</metric>
    <metric>upper_op2_start</metric>
    <metric>betweenness_final</metric>
    <metric>eigenvector_final</metric>
    <metric>clustering_final</metric>
    <metric>mean_path_final</metric>
    <metric>normalized_polarization_final</metric>
    <metric>ESBSG_polarization_final</metric>
    <metric>mean_op1_final</metric>
    <metric>mean_op2_final</metric>
    <metric>sd_op1_final</metric>
    <metric>sd_op2_final</metric>
    <metric>median_op1_final</metric>
    <metric>median_op2_final</metric>
    <metric>lower_op1_final</metric>
    <metric>lower_op2_final</metric>
    <metric>upper_op1_final</metric>
    <metric>upper_op2_final</metric>
    <steppedValueSet variable="RS" first="4" step="1" last="9"/>
    <enumeratedValueSet variable="cut-links-randomly?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="opinions">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="updating">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary">
      <value value="0.22"/>
      <value value="0.28"/>
      <value value="0.34"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="polar_repeats">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Y-opinion">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="draw_id_threshold">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="&quot;HK&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="3650"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Centroids_change">
      <value value="1.0E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-random">
      <value value="0.27"/>
      <value value="0.34"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tolerance-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-neis">
      <value value="64"/>
      <value value="52"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="use_opponents_ratio?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="polarisation-each-n-steps">
      <value value="10010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mode">
      <value value="&quot;openly-listen&quot;"/>
      <value value="&quot;vaguely-speak&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="d_threshold">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-each-n-steps">
      <value value="3650"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ESBG_furthest_out">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="identity_levels">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-length">
      <value value="40"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="id_threshold">
      <value value="0.39"/>
      <value value="0.49"/>
      <value value="0.59"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="network-change">
      <value value="&quot;link&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tolerance-level">
      <value value="1.08"/>
      <value value="0.648"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking-level">
      <value value="0.43"/>
      <value value="0.543"/>
      <value value="0.65"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="centroid_color?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="use_identity?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conformity-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N_centroids">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dissatisfied_updates_opinion">
      <value value="0.41"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="create-links-randomly?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="killing_centroids?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N-agents">
      <value value="129"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="construct-name?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="avoid-redundancies?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="HK-benchmark?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-comm-neis">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conformity-level">
      <value value="0.36"/>
      <value value="0.45"/>
      <value value="0.54"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold_drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="identity_type">
      <value value="&quot;global&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X-opinion">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="regularExperimentPart01LONG03" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="3651"/>
    <metric>betweenness_start</metric>
    <metric>eigenvector_start</metric>
    <metric>clustering_start</metric>
    <metric>mean_path_start</metric>
    <metric>normalized_polarization_start</metric>
    <metric>ESBSG_polarization_start</metric>
    <metric>mean_op1_start</metric>
    <metric>mean_op2_start</metric>
    <metric>sd_op1_start</metric>
    <metric>sd_op2_start</metric>
    <metric>median_op1_start</metric>
    <metric>median_op2_start</metric>
    <metric>lower_op1_start</metric>
    <metric>lower_op2_start</metric>
    <metric>upper_op1_start</metric>
    <metric>upper_op2_start</metric>
    <metric>betweenness_final</metric>
    <metric>eigenvector_final</metric>
    <metric>clustering_final</metric>
    <metric>mean_path_final</metric>
    <metric>normalized_polarization_final</metric>
    <metric>ESBSG_polarization_final</metric>
    <metric>mean_op1_final</metric>
    <metric>mean_op2_final</metric>
    <metric>sd_op1_final</metric>
    <metric>sd_op2_final</metric>
    <metric>median_op1_final</metric>
    <metric>median_op2_final</metric>
    <metric>lower_op1_final</metric>
    <metric>lower_op2_final</metric>
    <metric>upper_op1_final</metric>
    <metric>upper_op2_final</metric>
    <steppedValueSet variable="RS" first="10" step="1" last="19"/>
    <enumeratedValueSet variable="cut-links-randomly?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="opinions">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="updating">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary">
      <value value="0.22"/>
      <value value="0.28"/>
      <value value="0.34"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="polar_repeats">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Y-opinion">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="draw_id_threshold">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="&quot;HK&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="3650"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Centroids_change">
      <value value="1.0E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-random">
      <value value="0.27"/>
      <value value="0.34"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tolerance-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-neis">
      <value value="64"/>
      <value value="52"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="use_opponents_ratio?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="polarisation-each-n-steps">
      <value value="10010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mode">
      <value value="&quot;openly-listen&quot;"/>
      <value value="&quot;vaguely-speak&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="d_threshold">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-each-n-steps">
      <value value="3650"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ESBG_furthest_out">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="identity_levels">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-length">
      <value value="40"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="id_threshold">
      <value value="0.39"/>
      <value value="0.49"/>
      <value value="0.59"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="network-change">
      <value value="&quot;link&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tolerance-level">
      <value value="1.08"/>
      <value value="0.648"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking-level">
      <value value="0.43"/>
      <value value="0.543"/>
      <value value="0.65"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="centroid_color?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="use_identity?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conformity-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N_centroids">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dissatisfied_updates_opinion">
      <value value="0.41"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="create-links-randomly?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="killing_centroids?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N-agents">
      <value value="129"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="construct-name?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="avoid-redundancies?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="HK-benchmark?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-comm-neis">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conformity-level">
      <value value="0.36"/>
      <value value="0.45"/>
      <value value="0.54"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold_drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="identity_type">
      <value value="&quot;global&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X-opinion">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="regularExperimentPart01LONG04" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="3651"/>
    <metric>betweenness_start</metric>
    <metric>eigenvector_start</metric>
    <metric>clustering_start</metric>
    <metric>mean_path_start</metric>
    <metric>normalized_polarization_start</metric>
    <metric>ESBSG_polarization_start</metric>
    <metric>mean_op1_start</metric>
    <metric>mean_op2_start</metric>
    <metric>sd_op1_start</metric>
    <metric>sd_op2_start</metric>
    <metric>median_op1_start</metric>
    <metric>median_op2_start</metric>
    <metric>lower_op1_start</metric>
    <metric>lower_op2_start</metric>
    <metric>upper_op1_start</metric>
    <metric>upper_op2_start</metric>
    <metric>betweenness_final</metric>
    <metric>eigenvector_final</metric>
    <metric>clustering_final</metric>
    <metric>mean_path_final</metric>
    <metric>normalized_polarization_final</metric>
    <metric>ESBSG_polarization_final</metric>
    <metric>mean_op1_final</metric>
    <metric>mean_op2_final</metric>
    <metric>sd_op1_final</metric>
    <metric>sd_op2_final</metric>
    <metric>median_op1_final</metric>
    <metric>median_op2_final</metric>
    <metric>lower_op1_final</metric>
    <metric>lower_op2_final</metric>
    <metric>upper_op1_final</metric>
    <metric>upper_op2_final</metric>
    <steppedValueSet variable="RS" first="20" step="1" last="30"/>
    <enumeratedValueSet variable="cut-links-randomly?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="opinions">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="updating">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary">
      <value value="0.22"/>
      <value value="0.28"/>
      <value value="0.34"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="polar_repeats">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Y-opinion">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="draw_id_threshold">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="&quot;HK&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="3650"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Centroids_change">
      <value value="1.0E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-random">
      <value value="0.27"/>
      <value value="0.34"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tolerance-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-neis">
      <value value="64"/>
      <value value="52"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="use_opponents_ratio?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="polarisation-each-n-steps">
      <value value="10010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mode">
      <value value="&quot;openly-listen&quot;"/>
      <value value="&quot;vaguely-speak&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="d_threshold">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-each-n-steps">
      <value value="3650"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ESBG_furthest_out">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="identity_levels">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-length">
      <value value="40"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="id_threshold">
      <value value="0.39"/>
      <value value="0.49"/>
      <value value="0.59"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="network-change">
      <value value="&quot;link&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tolerance-level">
      <value value="1.08"/>
      <value value="0.648"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking-level">
      <value value="0.43"/>
      <value value="0.543"/>
      <value value="0.65"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="centroid_color?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="use_identity?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conformity-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N_centroids">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dissatisfied_updates_opinion">
      <value value="0.41"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="create-links-randomly?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="killing_centroids?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N-agents">
      <value value="129"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="construct-name?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="avoid-redundancies?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="HK-benchmark?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-comm-neis">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conformity-level">
      <value value="0.36"/>
      <value value="0.45"/>
      <value value="0.54"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold_drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="identity_type">
      <value value="&quot;global&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X-opinion">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="regularExperimentPart11" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="366"/>
    <metric>betweenness_start</metric>
    <metric>eigenvector_start</metric>
    <metric>clustering_start</metric>
    <metric>mean_path_start</metric>
    <metric>normalized_polarization_start</metric>
    <metric>ESBSG_polarization_start</metric>
    <metric>mean_op1_start</metric>
    <metric>mean_op2_start</metric>
    <metric>sd_op1_start</metric>
    <metric>sd_op2_start</metric>
    <metric>median_op1_start</metric>
    <metric>median_op2_start</metric>
    <metric>lower_op1_start</metric>
    <metric>lower_op2_start</metric>
    <metric>upper_op1_start</metric>
    <metric>upper_op2_start</metric>
    <metric>betweenness_final</metric>
    <metric>eigenvector_final</metric>
    <metric>clustering_final</metric>
    <metric>mean_path_final</metric>
    <metric>normalized_polarization_final</metric>
    <metric>ESBSG_polarization_final</metric>
    <metric>mean_op1_final</metric>
    <metric>mean_op2_final</metric>
    <metric>sd_op1_final</metric>
    <metric>sd_op2_final</metric>
    <metric>median_op1_final</metric>
    <metric>median_op2_final</metric>
    <metric>lower_op1_final</metric>
    <metric>lower_op2_final</metric>
    <metric>upper_op1_final</metric>
    <metric>upper_op2_final</metric>
    <steppedValueSet variable="RS" first="61" step="1" last="150"/>
    <enumeratedValueSet variable="cut-links-randomly?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="opinions">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="updating">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary">
      <value value="0.22"/>
      <value value="0.28"/>
      <value value="0.34"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="polar_repeats">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Y-opinion">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="draw_id_threshold">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="&quot;HK&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="365"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Centroids_change">
      <value value="1.0E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-random">
      <value value="0.27"/>
      <value value="0.34"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tolerance-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-neis">
      <value value="64"/>
      <value value="52"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="use_opponents_ratio?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="polarisation-each-n-steps">
      <value value="1010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mode">
      <value value="&quot;openly-listen&quot;"/>
      <value value="&quot;vaguely-speak&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="d_threshold">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-each-n-steps">
      <value value="365"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ESBG_furthest_out">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="identity_levels">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-length">
      <value value="40"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="id_threshold">
      <value value="0.39"/>
      <value value="0.49"/>
      <value value="0.59"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="network-change">
      <value value="&quot;link&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tolerance-level">
      <value value="1.08"/>
      <value value="0.648"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking-level">
      <value value="0.43"/>
      <value value="0.543"/>
      <value value="0.65"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="centroid_color?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="use_identity?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conformity-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N_centroids">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dissatisfied_updates_opinion">
      <value value="0.41"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="create-links-randomly?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="killing_centroids?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N-agents">
      <value value="129"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="construct-name?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="avoid-redundancies?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="HK-benchmark?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-comm-neis">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conformity-level">
      <value value="0.36"/>
      <value value="0.45"/>
      <value value="0.54"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold_drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="identity_type">
      <value value="&quot;global&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X-opinion">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="regularExperimentPart13" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="366"/>
    <metric>betweenness_start</metric>
    <metric>eigenvector_start</metric>
    <metric>clustering_start</metric>
    <metric>mean_path_start</metric>
    <metric>normalized_polarization_start</metric>
    <metric>ESBSG_polarization_start</metric>
    <metric>mean_op1_start</metric>
    <metric>mean_op2_start</metric>
    <metric>sd_op1_start</metric>
    <metric>sd_op2_start</metric>
    <metric>median_op1_start</metric>
    <metric>median_op2_start</metric>
    <metric>lower_op1_start</metric>
    <metric>lower_op2_start</metric>
    <metric>upper_op1_start</metric>
    <metric>upper_op2_start</metric>
    <metric>betweenness_final</metric>
    <metric>eigenvector_final</metric>
    <metric>clustering_final</metric>
    <metric>mean_path_final</metric>
    <metric>normalized_polarization_final</metric>
    <metric>ESBSG_polarization_final</metric>
    <metric>mean_op1_final</metric>
    <metric>mean_op2_final</metric>
    <metric>sd_op1_final</metric>
    <metric>sd_op2_final</metric>
    <metric>median_op1_final</metric>
    <metric>median_op2_final</metric>
    <metric>lower_op1_final</metric>
    <metric>lower_op2_final</metric>
    <metric>upper_op1_final</metric>
    <metric>upper_op2_final</metric>
    <steppedValueSet variable="RS" first="61" step="1" last="150"/>
    <enumeratedValueSet variable="cut-links-randomly?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="opinions">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="updating">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary">
      <value value="0.22"/>
      <value value="0.28"/>
      <value value="0.34"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="polar_repeats">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Y-opinion">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="draw_id_threshold">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="&quot;HK&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="365"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Centroids_change">
      <value value="1.0E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-random">
      <value value="0.27"/>
      <value value="0.34"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tolerance-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-neis">
      <value value="64"/>
      <value value="52"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="use_opponents_ratio?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="polarisation-each-n-steps">
      <value value="1010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mode">
      <value value="&quot;openly-listen&quot;"/>
      <value value="&quot;vaguely-speak&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="d_threshold">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-each-n-steps">
      <value value="365"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ESBG_furthest_out">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="identity_levels">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-length">
      <value value="40"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="id_threshold">
      <value value="0.49"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="network-change">
      <value value="&quot;link&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tolerance-level">
      <value value="1.08"/>
      <value value="0.648"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking-level">
      <value value="0.43"/>
      <value value="0.543"/>
      <value value="0.65"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="centroid_color?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="use_identity?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conformity-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N_centroids">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dissatisfied_updates_opinion">
      <value value="0.41"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="create-links-randomly?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="killing_centroids?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N-agents">
      <value value="129"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="construct-name?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="avoid-redundancies?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="HK-benchmark?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-comm-neis">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conformity-level">
      <value value="0.36"/>
      <value value="0.45"/>
      <value value="0.54"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold_drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="identity_type">
      <value value="&quot;global&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X-opinion">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="regularExperimentPart14" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="366"/>
    <metric>betweenness_start</metric>
    <metric>eigenvector_start</metric>
    <metric>clustering_start</metric>
    <metric>mean_path_start</metric>
    <metric>normalized_polarization_start</metric>
    <metric>ESBSG_polarization_start</metric>
    <metric>mean_op1_start</metric>
    <metric>mean_op2_start</metric>
    <metric>sd_op1_start</metric>
    <metric>sd_op2_start</metric>
    <metric>median_op1_start</metric>
    <metric>median_op2_start</metric>
    <metric>lower_op1_start</metric>
    <metric>lower_op2_start</metric>
    <metric>upper_op1_start</metric>
    <metric>upper_op2_start</metric>
    <metric>betweenness_final</metric>
    <metric>eigenvector_final</metric>
    <metric>clustering_final</metric>
    <metric>mean_path_final</metric>
    <metric>normalized_polarization_final</metric>
    <metric>ESBSG_polarization_final</metric>
    <metric>mean_op1_final</metric>
    <metric>mean_op2_final</metric>
    <metric>sd_op1_final</metric>
    <metric>sd_op2_final</metric>
    <metric>median_op1_final</metric>
    <metric>median_op2_final</metric>
    <metric>lower_op1_final</metric>
    <metric>lower_op2_final</metric>
    <metric>upper_op1_final</metric>
    <metric>upper_op2_final</metric>
    <steppedValueSet variable="RS" first="1" step="1" last="60"/>
    <enumeratedValueSet variable="n-neis">
      <value value="52"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="use_identity?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tolerance-level">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cut-links-randomly?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="opinions">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="updating">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary">
      <value value="0.22"/>
      <value value="0.28"/>
      <value value="0.34"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="polar_repeats">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Y-opinion">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="draw_id_threshold">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="&quot;HK&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="365"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Centroids_change">
      <value value="1.0E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-random">
      <value value="0.27"/>
      <value value="0.34"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tolerance-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="use_opponents_ratio?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="polarisation-each-n-steps">
      <value value="1010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mode">
      <value value="&quot;openly-listen&quot;"/>
      <value value="&quot;vaguely-speak&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="d_threshold">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-each-n-steps">
      <value value="365"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ESBG_furthest_out">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="identity_levels">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-length">
      <value value="40"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="id_threshold">
      <value value="0.49"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="network-change">
      <value value="&quot;link&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking-level">
      <value value="0.43"/>
      <value value="0.543"/>
      <value value="0.65"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="centroid_color?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conformity-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N_centroids">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dissatisfied_updates_opinion">
      <value value="0.41"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="create-links-randomly?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="killing_centroids?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N-agents">
      <value value="129"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="construct-name?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="avoid-redundancies?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="HK-benchmark?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-comm-neis">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conformity-level">
      <value value="0.36"/>
      <value value="0.45"/>
      <value value="0.54"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold_drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="identity_type">
      <value value="&quot;global&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X-opinion">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="regularExperimentPart15" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="366"/>
    <metric>betweenness_start</metric>
    <metric>eigenvector_start</metric>
    <metric>clustering_start</metric>
    <metric>mean_path_start</metric>
    <metric>normalized_polarization_start</metric>
    <metric>ESBSG_polarization_start</metric>
    <metric>mean_op1_start</metric>
    <metric>mean_op2_start</metric>
    <metric>sd_op1_start</metric>
    <metric>sd_op2_start</metric>
    <metric>median_op1_start</metric>
    <metric>median_op2_start</metric>
    <metric>lower_op1_start</metric>
    <metric>lower_op2_start</metric>
    <metric>upper_op1_start</metric>
    <metric>upper_op2_start</metric>
    <metric>betweenness_final</metric>
    <metric>eigenvector_final</metric>
    <metric>clustering_final</metric>
    <metric>mean_path_final</metric>
    <metric>normalized_polarization_final</metric>
    <metric>ESBSG_polarization_final</metric>
    <metric>mean_op1_final</metric>
    <metric>mean_op2_final</metric>
    <metric>sd_op1_final</metric>
    <metric>sd_op2_final</metric>
    <metric>median_op1_final</metric>
    <metric>median_op2_final</metric>
    <metric>lower_op1_final</metric>
    <metric>lower_op2_final</metric>
    <metric>upper_op1_final</metric>
    <metric>upper_op2_final</metric>
    <steppedValueSet variable="RS" first="1" step="1" last="60"/>
    <enumeratedValueSet variable="n-neis">
      <value value="13"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="use_identity?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tolerance-level">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cut-links-randomly?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="opinions">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="updating">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary">
      <value value="0.22"/>
      <value value="0.28"/>
      <value value="0.34"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="polar_repeats">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Y-opinion">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="draw_id_threshold">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="&quot;HK&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="365"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Centroids_change">
      <value value="1.0E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-random">
      <value value="0.27"/>
      <value value="0.34"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tolerance-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="use_opponents_ratio?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="polarisation-each-n-steps">
      <value value="1010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mode">
      <value value="&quot;openly-listen&quot;"/>
      <value value="&quot;vaguely-speak&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="d_threshold">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-each-n-steps">
      <value value="365"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ESBG_furthest_out">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="identity_levels">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-length">
      <value value="40"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="id_threshold">
      <value value="0.49"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="network-change">
      <value value="&quot;link&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking-level">
      <value value="0.43"/>
      <value value="0.543"/>
      <value value="0.65"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="centroid_color?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conformity-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N_centroids">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dissatisfied_updates_opinion">
      <value value="0.41"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="create-links-randomly?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="killing_centroids?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N-agents">
      <value value="129"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="construct-name?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="avoid-redundancies?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="HK-benchmark?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-comm-neis">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conformity-level">
      <value value="0.36"/>
      <value value="0.45"/>
      <value value="0.54"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold_drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="identity_type">
      <value value="&quot;global&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X-opinion">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="regularExperimentPart16" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="366"/>
    <metric>betweenness_start</metric>
    <metric>eigenvector_start</metric>
    <metric>clustering_start</metric>
    <metric>mean_path_start</metric>
    <metric>normalized_polarization_start</metric>
    <metric>ESBSG_polarization_start</metric>
    <metric>mean_op1_start</metric>
    <metric>mean_op2_start</metric>
    <metric>sd_op1_start</metric>
    <metric>sd_op2_start</metric>
    <metric>median_op1_start</metric>
    <metric>median_op2_start</metric>
    <metric>lower_op1_start</metric>
    <metric>lower_op2_start</metric>
    <metric>upper_op1_start</metric>
    <metric>upper_op2_start</metric>
    <metric>betweenness_final</metric>
    <metric>eigenvector_final</metric>
    <metric>clustering_final</metric>
    <metric>mean_path_final</metric>
    <metric>normalized_polarization_final</metric>
    <metric>ESBSG_polarization_final</metric>
    <metric>mean_op1_final</metric>
    <metric>mean_op2_final</metric>
    <metric>sd_op1_final</metric>
    <metric>sd_op2_final</metric>
    <metric>median_op1_final</metric>
    <metric>median_op2_final</metric>
    <metric>lower_op1_final</metric>
    <metric>lower_op2_final</metric>
    <metric>upper_op1_final</metric>
    <metric>upper_op2_final</metric>
    <steppedValueSet variable="RS" first="1" step="1" last="60"/>
    <enumeratedValueSet variable="n-neis">
      <value value="13"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="use_identity?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tolerance-level">
      <value value="0.648"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cut-links-randomly?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="opinions">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="updating">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary">
      <value value="0.22"/>
      <value value="0.28"/>
      <value value="0.34"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="polar_repeats">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Y-opinion">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="draw_id_threshold">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="&quot;HK&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="365"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Centroids_change">
      <value value="1.0E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-random">
      <value value="0.27"/>
      <value value="0.34"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tolerance-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="use_opponents_ratio?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="polarisation-each-n-steps">
      <value value="1010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mode">
      <value value="&quot;openly-listen&quot;"/>
      <value value="&quot;vaguely-speak&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="d_threshold">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-each-n-steps">
      <value value="365"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ESBG_furthest_out">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="identity_levels">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-length">
      <value value="40"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="id_threshold">
      <value value="0.49"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="network-change">
      <value value="&quot;link&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking-level">
      <value value="0.43"/>
      <value value="0.543"/>
      <value value="0.65"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="centroid_color?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conformity-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N_centroids">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dissatisfied_updates_opinion">
      <value value="0.41"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="create-links-randomly?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="killing_centroids?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N-agents">
      <value value="129"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="construct-name?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="avoid-redundancies?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="HK-benchmark?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-comm-neis">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conformity-level">
      <value value="0.36"/>
      <value value="0.45"/>
      <value value="0.54"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold_drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="identity_type">
      <value value="&quot;global&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X-opinion">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="regularExperimentPart17a" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="366"/>
    <metric>betweenness_start</metric>
    <metric>eigenvector_start</metric>
    <metric>clustering_start</metric>
    <metric>mean_path_start</metric>
    <metric>normalized_polarization_start</metric>
    <metric>ESBSG_polarization_start</metric>
    <metric>mean_op1_start</metric>
    <metric>mean_op2_start</metric>
    <metric>sd_op1_start</metric>
    <metric>sd_op2_start</metric>
    <metric>median_op1_start</metric>
    <metric>median_op2_start</metric>
    <metric>lower_op1_start</metric>
    <metric>lower_op2_start</metric>
    <metric>upper_op1_start</metric>
    <metric>upper_op2_start</metric>
    <metric>betweenness_final</metric>
    <metric>eigenvector_final</metric>
    <metric>clustering_final</metric>
    <metric>mean_path_final</metric>
    <metric>normalized_polarization_final</metric>
    <metric>ESBSG_polarization_final</metric>
    <metric>mean_op1_final</metric>
    <metric>mean_op2_final</metric>
    <metric>sd_op1_final</metric>
    <metric>sd_op2_final</metric>
    <metric>median_op1_final</metric>
    <metric>median_op2_final</metric>
    <metric>lower_op1_final</metric>
    <metric>lower_op2_final</metric>
    <metric>upper_op1_final</metric>
    <metric>upper_op2_final</metric>
    <steppedValueSet variable="RS" first="1" step="1" last="150"/>
    <enumeratedValueSet variable="cut-links-randomly?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="opinions">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="updating">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary">
      <value value="0.22"/>
      <value value="0.28"/>
      <value value="0.34"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="polar_repeats">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Y-opinion">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="draw_id_threshold">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="&quot;HK&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="365"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Centroids_change">
      <value value="1.0E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-random">
      <value value="0.27"/>
      <value value="0.34"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tolerance-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-neis">
      <value value="64"/>
      <value value="52"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="use_opponents_ratio?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="polarisation-each-n-steps">
      <value value="1010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mode">
      <value value="&quot;openly-listen&quot;"/>
      <value value="&quot;vaguely-speak&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="d_threshold">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-each-n-steps">
      <value value="365"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ESBG_furthest_out">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="identity_levels">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-length">
      <value value="40"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="id_threshold">
      <value value="0.39"/>
      <value value="0.49"/>
      <value value="0.59"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="network-change">
      <value value="&quot;link&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tolerance-level">
      <value value="1.08"/>
      <value value="0.648"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking-level">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="centroid_color?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="use_identity?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conformity-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N_centroids">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dissatisfied_updates_opinion">
      <value value="0.41"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="create-links-randomly?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="killing_centroids?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N-agents">
      <value value="129"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="construct-name?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="avoid-redundancies?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="HK-benchmark?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-comm-neis">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conformity-level">
      <value value="0.36"/>
      <value value="0.45"/>
      <value value="0.54"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold_drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="identity_type">
      <value value="&quot;global&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X-opinion">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="regularExperimentPart17b" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="366"/>
    <metric>betweenness_start</metric>
    <metric>eigenvector_start</metric>
    <metric>clustering_start</metric>
    <metric>mean_path_start</metric>
    <metric>normalized_polarization_start</metric>
    <metric>ESBSG_polarization_start</metric>
    <metric>mean_op1_start</metric>
    <metric>mean_op2_start</metric>
    <metric>sd_op1_start</metric>
    <metric>sd_op2_start</metric>
    <metric>median_op1_start</metric>
    <metric>median_op2_start</metric>
    <metric>lower_op1_start</metric>
    <metric>lower_op2_start</metric>
    <metric>upper_op1_start</metric>
    <metric>upper_op2_start</metric>
    <metric>betweenness_final</metric>
    <metric>eigenvector_final</metric>
    <metric>clustering_final</metric>
    <metric>mean_path_final</metric>
    <metric>normalized_polarization_final</metric>
    <metric>ESBSG_polarization_final</metric>
    <metric>mean_op1_final</metric>
    <metric>mean_op2_final</metric>
    <metric>sd_op1_final</metric>
    <metric>sd_op2_final</metric>
    <metric>median_op1_final</metric>
    <metric>median_op2_final</metric>
    <metric>lower_op1_final</metric>
    <metric>lower_op2_final</metric>
    <metric>upper_op1_final</metric>
    <metric>upper_op2_final</metric>
    <steppedValueSet variable="RS" first="1" step="1" last="150"/>
    <enumeratedValueSet variable="cut-links-randomly?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="opinions">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="updating">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="polar_repeats">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Y-opinion">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="draw_id_threshold">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="&quot;HK&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="365"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Centroids_change">
      <value value="1.0E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-random">
      <value value="0.27"/>
      <value value="0.34"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tolerance-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-neis">
      <value value="64"/>
      <value value="52"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="use_opponents_ratio?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="polarisation-each-n-steps">
      <value value="1010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mode">
      <value value="&quot;openly-listen&quot;"/>
      <value value="&quot;vaguely-speak&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="d_threshold">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-each-n-steps">
      <value value="365"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ESBG_furthest_out">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="identity_levels">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-length">
      <value value="40"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="id_threshold">
      <value value="0.39"/>
      <value value="0.49"/>
      <value value="0.59"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="network-change">
      <value value="&quot;link&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tolerance-level">
      <value value="1.08"/>
      <value value="0.648"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking-level">
      <value value="0.43"/>
      <value value="0.543"/>
      <value value="0.65"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="centroid_color?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="use_identity?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conformity-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N_centroids">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dissatisfied_updates_opinion">
      <value value="0.41"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="create-links-randomly?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="killing_centroids?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N-agents">
      <value value="129"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="construct-name?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="avoid-redundancies?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="HK-benchmark?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-comm-neis">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conformity-level">
      <value value="0.36"/>
      <value value="0.45"/>
      <value value="0.54"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold_drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="identity_type">
      <value value="&quot;global&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X-opinion">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="regularExperimentPart18a" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="366"/>
    <metric>betweenness_start</metric>
    <metric>eigenvector_start</metric>
    <metric>clustering_start</metric>
    <metric>mean_path_start</metric>
    <metric>normalized_polarization_start</metric>
    <metric>ESBSG_polarization_start</metric>
    <metric>mean_op1_start</metric>
    <metric>mean_op2_start</metric>
    <metric>sd_op1_start</metric>
    <metric>sd_op2_start</metric>
    <metric>median_op1_start</metric>
    <metric>median_op2_start</metric>
    <metric>lower_op1_start</metric>
    <metric>lower_op2_start</metric>
    <metric>upper_op1_start</metric>
    <metric>upper_op2_start</metric>
    <metric>betweenness_final</metric>
    <metric>eigenvector_final</metric>
    <metric>clustering_final</metric>
    <metric>mean_path_final</metric>
    <metric>normalized_polarization_final</metric>
    <metric>ESBSG_polarization_final</metric>
    <metric>mean_op1_final</metric>
    <metric>mean_op2_final</metric>
    <metric>sd_op1_final</metric>
    <metric>sd_op2_final</metric>
    <metric>median_op1_final</metric>
    <metric>median_op2_final</metric>
    <metric>lower_op1_final</metric>
    <metric>lower_op2_final</metric>
    <metric>upper_op1_final</metric>
    <metric>upper_op2_final</metric>
    <steppedValueSet variable="RS" first="1" step="1" last="150"/>
    <enumeratedValueSet variable="cut-links-randomly?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="opinions">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="updating">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary">
      <value value="0.22"/>
      <value value="0.28"/>
      <value value="0.34"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="polar_repeats">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Y-opinion">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="draw_id_threshold">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="&quot;HK&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="365"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Centroids_change">
      <value value="1.0E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-random">
      <value value="0.27"/>
      <value value="0.34"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tolerance-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-neis">
      <value value="64"/>
      <value value="52"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="use_opponents_ratio?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="polarisation-each-n-steps">
      <value value="1010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mode">
      <value value="&quot;openly-listen&quot;"/>
      <value value="&quot;vaguely-speak&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="d_threshold">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-each-n-steps">
      <value value="365"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ESBG_furthest_out">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="identity_levels">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-length">
      <value value="40"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="id_threshold">
      <value value="0.49"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="network-change">
      <value value="&quot;link&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tolerance-level">
      <value value="1.08"/>
      <value value="0.648"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking-level">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="centroid_color?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="use_identity?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conformity-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N_centroids">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dissatisfied_updates_opinion">
      <value value="0.41"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="create-links-randomly?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="killing_centroids?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N-agents">
      <value value="129"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="construct-name?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="avoid-redundancies?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="HK-benchmark?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-comm-neis">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conformity-level">
      <value value="0.36"/>
      <value value="0.45"/>
      <value value="0.54"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold_drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="identity_type">
      <value value="&quot;global&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X-opinion">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="regularExperimentPart18b" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="366"/>
    <metric>betweenness_start</metric>
    <metric>eigenvector_start</metric>
    <metric>clustering_start</metric>
    <metric>mean_path_start</metric>
    <metric>normalized_polarization_start</metric>
    <metric>ESBSG_polarization_start</metric>
    <metric>mean_op1_start</metric>
    <metric>mean_op2_start</metric>
    <metric>sd_op1_start</metric>
    <metric>sd_op2_start</metric>
    <metric>median_op1_start</metric>
    <metric>median_op2_start</metric>
    <metric>lower_op1_start</metric>
    <metric>lower_op2_start</metric>
    <metric>upper_op1_start</metric>
    <metric>upper_op2_start</metric>
    <metric>betweenness_final</metric>
    <metric>eigenvector_final</metric>
    <metric>clustering_final</metric>
    <metric>mean_path_final</metric>
    <metric>normalized_polarization_final</metric>
    <metric>ESBSG_polarization_final</metric>
    <metric>mean_op1_final</metric>
    <metric>mean_op2_final</metric>
    <metric>sd_op1_final</metric>
    <metric>sd_op2_final</metric>
    <metric>median_op1_final</metric>
    <metric>median_op2_final</metric>
    <metric>lower_op1_final</metric>
    <metric>lower_op2_final</metric>
    <metric>upper_op1_final</metric>
    <metric>upper_op2_final</metric>
    <steppedValueSet variable="RS" first="1" step="1" last="150"/>
    <enumeratedValueSet variable="cut-links-randomly?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="opinions">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="updating">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="polar_repeats">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Y-opinion">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="draw_id_threshold">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="&quot;HK&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="365"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Centroids_change">
      <value value="1.0E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-random">
      <value value="0.27"/>
      <value value="0.34"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tolerance-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-neis">
      <value value="64"/>
      <value value="52"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="use_opponents_ratio?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="polarisation-each-n-steps">
      <value value="1010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mode">
      <value value="&quot;openly-listen&quot;"/>
      <value value="&quot;vaguely-speak&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="d_threshold">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-each-n-steps">
      <value value="365"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ESBG_furthest_out">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="identity_levels">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-length">
      <value value="40"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="id_threshold">
      <value value="0.49"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="network-change">
      <value value="&quot;link&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tolerance-level">
      <value value="1.08"/>
      <value value="0.648"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking-level">
      <value value="0.43"/>
      <value value="0.543"/>
      <value value="0.65"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="centroid_color?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="use_identity?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conformity-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N_centroids">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dissatisfied_updates_opinion">
      <value value="0.41"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="create-links-randomly?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="killing_centroids?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N-agents">
      <value value="129"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="construct-name?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="avoid-redundancies?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="HK-benchmark?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-comm-neis">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conformity-level">
      <value value="0.36"/>
      <value value="0.45"/>
      <value value="0.54"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold_drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="identity_type">
      <value value="&quot;global&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X-opinion">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="regularExperimentPart31" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="366"/>
    <metric>betweenness_start</metric>
    <metric>eigenvector_start</metric>
    <metric>clustering_start</metric>
    <metric>mean_path_start</metric>
    <metric>normalized_polarization_start</metric>
    <metric>ESBSG_polarization_start</metric>
    <metric>mean_op1_start</metric>
    <metric>mean_op2_start</metric>
    <metric>mean_op3_start</metric>
    <metric>mean_op4_start</metric>
    <metric>sd_op1_start</metric>
    <metric>sd_op2_start</metric>
    <metric>sd_op3_start</metric>
    <metric>sd_op4_start</metric>
    <metric>median_op1_start</metric>
    <metric>median_op2_start</metric>
    <metric>median_op3_start</metric>
    <metric>median_op4_start</metric>
    <metric>lower_op1_start</metric>
    <metric>lower_op2_start</metric>
    <metric>lower_op3_start</metric>
    <metric>lower_op4_start</metric>
    <metric>upper_op1_start</metric>
    <metric>upper_op2_start</metric>
    <metric>upper_op3_start</metric>
    <metric>upper_op4_start</metric>
    <metric>betweenness_final</metric>
    <metric>eigenvector_final</metric>
    <metric>clustering_final</metric>
    <metric>mean_path_final</metric>
    <metric>normalized_polarization_final</metric>
    <metric>ESBSG_polarization_final</metric>
    <metric>mean_op1_final</metric>
    <metric>mean_op2_final</metric>
    <metric>mean_op3_final</metric>
    <metric>mean_op4_final</metric>
    <metric>sd_op1_final</metric>
    <metric>sd_op2_final</metric>
    <metric>sd_op3_final</metric>
    <metric>sd_op4_final</metric>
    <metric>median_op1_final</metric>
    <metric>median_op2_final</metric>
    <metric>median_op3_final</metric>
    <metric>median_op4_final</metric>
    <metric>lower_op1_final</metric>
    <metric>lower_op2_final</metric>
    <metric>lower_op3_final</metric>
    <metric>lower_op4_final</metric>
    <metric>upper_op1_final</metric>
    <metric>upper_op2_final</metric>
    <metric>upper_op3_final</metric>
    <metric>upper_op4_final</metric>
    <steppedValueSet variable="RS" first="1" step="1" last="150"/>
    <enumeratedValueSet variable="use_identity?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="opinions">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary">
      <value value="0.1"/>
      <value value="0.22"/>
      <value value="0.28"/>
      <value value="0.34"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conformity-level">
      <value value="0.45"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="id_threshold">
      <value value="0.39"/>
      <value value="0.49"/>
      <value value="0.59"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tolerance-level">
      <value value="1.08"/>
      <value value="0.648"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking-level">
      <value value="0.5"/>
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mode">
      <value value="&quot;openly-listen&quot;"/>
      <value value="&quot;vaguely-speak&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-random">
      <value value="0.27"/>
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-neis">
      <value value="52"/>
      <value value="13"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="365"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-each-n-steps">
      <value value="365"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="polar_repeats">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="draw_id_threshold">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="&quot;HK&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Centroids_change">
      <value value="1.0E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tolerance-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="use_opponents_ratio?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="polarisation-each-n-steps">
      <value value="1010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="d_threshold">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ESBG_furthest_out">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="identity_levels">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-length">
      <value value="40"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="centroid_color?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conformity-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N_centroids">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dissatisfied_updates_opinion">
      <value value="0.41"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="create-links-randomly?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="killing_centroids?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N-agents">
      <value value="129"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="construct-name?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="avoid-redundancies?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="HK-benchmark?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-comm-neis">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold_drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="identity_type">
      <value value="&quot;global&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X-opinion">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Y-opinion">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cut-links-randomly?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="updating">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="network-change">
      <value value="&quot;link&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="regularExperimentPart33" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="366"/>
    <metric>betweenness_start</metric>
    <metric>eigenvector_start</metric>
    <metric>clustering_start</metric>
    <metric>mean_path_start</metric>
    <metric>normalized_polarization_start</metric>
    <metric>ESBSG_polarization_start</metric>
    <metric>mean_op1_start</metric>
    <metric>mean_op2_start</metric>
    <metric>mean_op3_start</metric>
    <metric>mean_op4_start</metric>
    <metric>sd_op1_start</metric>
    <metric>sd_op2_start</metric>
    <metric>sd_op3_start</metric>
    <metric>sd_op4_start</metric>
    <metric>median_op1_start</metric>
    <metric>median_op2_start</metric>
    <metric>median_op3_start</metric>
    <metric>median_op4_start</metric>
    <metric>lower_op1_start</metric>
    <metric>lower_op2_start</metric>
    <metric>lower_op3_start</metric>
    <metric>lower_op4_start</metric>
    <metric>upper_op1_start</metric>
    <metric>upper_op2_start</metric>
    <metric>upper_op3_start</metric>
    <metric>upper_op4_start</metric>
    <metric>betweenness_final</metric>
    <metric>eigenvector_final</metric>
    <metric>clustering_final</metric>
    <metric>mean_path_final</metric>
    <metric>normalized_polarization_final</metric>
    <metric>ESBSG_polarization_final</metric>
    <metric>mean_op1_final</metric>
    <metric>mean_op2_final</metric>
    <metric>mean_op3_final</metric>
    <metric>mean_op4_final</metric>
    <metric>sd_op1_final</metric>
    <metric>sd_op2_final</metric>
    <metric>sd_op3_final</metric>
    <metric>sd_op4_final</metric>
    <metric>median_op1_final</metric>
    <metric>median_op2_final</metric>
    <metric>median_op3_final</metric>
    <metric>median_op4_final</metric>
    <metric>lower_op1_final</metric>
    <metric>lower_op2_final</metric>
    <metric>lower_op3_final</metric>
    <metric>lower_op4_final</metric>
    <metric>upper_op1_final</metric>
    <metric>upper_op2_final</metric>
    <metric>upper_op3_final</metric>
    <metric>upper_op4_final</metric>
    <steppedValueSet variable="RS" first="1" step="1" last="150"/>
    <enumeratedValueSet variable="use_identity?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="opinions">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary">
      <value value="0.1"/>
      <value value="0.22"/>
      <value value="0.28"/>
      <value value="0.34"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conformity-level">
      <value value="0.45"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="id_threshold">
      <value value="0.49"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tolerance-level">
      <value value="1.08"/>
      <value value="0.648"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking-level">
      <value value="0.5"/>
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mode">
      <value value="&quot;openly-listen&quot;"/>
      <value value="&quot;vaguely-speak&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-random">
      <value value="0.27"/>
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-neis">
      <value value="52"/>
      <value value="13"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="365"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-each-n-steps">
      <value value="365"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="polar_repeats">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="draw_id_threshold">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="&quot;HK&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Centroids_change">
      <value value="1.0E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tolerance-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="use_opponents_ratio?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="polarisation-each-n-steps">
      <value value="1010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="d_threshold">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ESBG_furthest_out">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="identity_levels">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-length">
      <value value="40"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="centroid_color?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conformity-drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N_centroids">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dissatisfied_updates_opinion">
      <value value="0.41"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="create-links-randomly?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="killing_centroids?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N-agents">
      <value value="129"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="construct-name?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="avoid-redundancies?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="HK-benchmark?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-comm-neis">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold_drawn">
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="identity_type">
      <value value="&quot;global&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X-opinion">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Y-opinion">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cut-links-randomly?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="updating">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="network-change">
      <value value="&quot;link&quot;"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
