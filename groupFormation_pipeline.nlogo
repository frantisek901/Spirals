;; Model for study of mutual interaction between public opinion and network structure.
;; RQ: How polarisation of public opinion shapes the network structure
;;     and how the network structure shapes public opinion and influences individual opinions?
;;
;; This code is from project of Mike, Ashley, Ashwin and FranCesko,
;; we apply Hegselmann-Krausse model in more than 1D and look how agents adapt in >1D opinion space and whether they form groups,
;; then we include small-world network (Watts-Strogatz) as another constraint, features of spiral of silence, and
;; individual assignment of uncertainity, tollerance, comformity and outspokeness.
;;
;; !!!FROZEN FOR PIPELINE!!!
;;

;; Created:  2021-10-21 FranCesko
;; Edited:   2022-02-04 FranCesko
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

extensions [nw]

breed [centroids centroid]
undirected-link-breed [comms comm]
undirected-link-breed [l-distances l-distance]

turtles-own [Opinion-position P-speaking Speak? Uncertainty Record Last-opinion Pol-bias Initial-opinion Tolerance Conformity Satisfied? group distance_to_centroid]
l-distances-own [l-weight]

globals [main-Record components positions network-changes agents polarisation normalized_polarisation unweighted_polarisation unweighted_normalized_polarisation]


;; Initialization and setup
to setup
  ;; Redundant conditions which should be avoided -- if the boundary is drawn as constant, then is completely same whether agents vaguely speak or openly listen,
  ;; seme case is for probability speaking of 100%, then it's same whether individual probability is drawn as constant or uniform, result is still same: all agents has probability 100%.
  if avoid-redundancies? and mode = "vaguely-speak" and boundary-drawn = "constant" [stop]
  if avoid-redundancies? and p-speaking-level = 1 and p-speaking-drawn = "uniform" [stop]
  ;; these two conditions cover 7/16 of all simulations, approx. the half! This code should stop them from running.
  ;(avoid-redundancies? and mode = "vaguely-speak" and boundary-drawn = "constant") or (avoid-redundancies? and p-speaking-level = 1 and p-speaking-drawn = "uniform")

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
    set Conformity get-conformity  ;; setting individual conformity level.
    getColor  ;; Coloring the agents according their opinion.
    getPlace  ;; Moving agents to the opinion space according their opinions.
  ]
  set agents turtle-set turtles  ;; Note: If we just write 'set agents turtles', then variable 'agents' is a synonym for 'turtles', so it will contain in the future created centroids!


  ;; Coloring patches according the number of agents/turtles on them.
  ask patches [set pcolor patch-color]
  ;; Hiding links so to improve simulation speed performance.
  ask comms [set hidden? TRUE]
  ;; Setting the indicator of change for the whole simulation, again as non-stable.
  set main-Record n-values record-length [0]
  ;; Setting control variable of network changes
  set network-changes 0
  ;; Compute polarisation
  compute-polarisation-repeatedly

  reset-ticks

  ;;;; Finally, we record initial state of simulation
  ;; If we want we could construct filename to contain all important parameters shaping initial condition, so the name is unique stamp of initial state!
  if construct-name? [set file-name-core (word RS "_" N-agents "_" p-random "_" n-neis "_" opinions "_" updating "_" boundary "_" boundary-drawn "_" p-speaking-level "_"  p-speaking-drawn "_" mode)]
  ;; recording itself
  if record? [record-state-of-simulation]
end

;; We compute polarisation several times and then set it for the average
to compute-polarisation-repeatedly
  ;; Initialization of temporal variables
  let r 0
  let p []
  let np []
  let up []
  let unp []

  ;; Repeating cycle
  while [r < polar_repeats] [
    compute-polarisation
    set p lput polarisation p
    set np lput normalized_polarisation np
    set up lput unweighted_polarisation up
    set unp lput unweighted_normalized_polarisation unp
    set r r + 1
  ]

  ;; Setting variables back
  set polarisation precision (mean p) 3
  set normalized_polarisation precision (mean np) 3
  set unweighted_polarisation precision (mean up) 3
  set unweighted_normalized_polarisation precision (mean unp) 3
end

;; NOTE: Now I am iplementing it for N = 2 centroids, but I prepare code for easy generalisation for N > 2.
to compute-polarisation
  ;; Preparation
  ask centroids [die]
  ;let original_centroids_value N_centroids

  ;; Detection of clusters via Louvain
  ask l-distances [die]  ;; Cleaning environment
  ask agents [create-l-distances-with other agents with [(sqrt(4 * opinions) - opinion-distance) / sqrt(4 * opinions) >= id_threshold] [set l-weight opinion-distance2 ([opinion-position] of end1)([opinion-position] of end2)]]
  ;show count l-distances
  ;ask l-distances with [l-weight < d_threshold] [die]
  nw:set-context agents l-distances ;with [l-weight >= d_threshold]
  let communities nw:louvain-communities
  ;show count l-distances
  ask l-distances [die]
  set N_centroids length communities

  ;; Computing clusters' mean 'opinion-position'
  let postions-clusters [] ;; List with all positions of all clusters
  foreach communities [c ->
    let one []  ;; List for one positio nof one cluster
    foreach range opinions [o ->
      set one lput precision (mean [item o opinion-position] of c) 3 one
    ]
    ;show one
    set postions-clusters lput one postions-clusters
    ;show postions-clusters
  ]

  ;; Preparation of centroids -- feedeing them with communities
  create-centroids N_centroids [
    set heading (who - min [who] of centroids)
    set Opinion-position item heading postions-clusters  ;; We set opinions, we try to do it smoothly...
    ;show Opinion-position
    set shape "circle"
    set size 1.5
    set color 5 + who * 10
    getPlace
  ]

  ;; Assignment of agents to groups
  ;let min_group min [who] of centroids
  ask agents [set group [who] of min-one-of centroids [opinion-distance]]

  ;; Computation of centroids possitions
  compute-centroids-positions

  ;let iter 0

  ;; Iterating cycle -- looking for good match of centroids
  while [sum [opinion-distance3 (Last-opinion) (Opinion-position)] of centroids > Centroids_change] [

    ;set iter iter + 1
    ;show (word "Iteration: " iter)

    ;; turtles compute whether they are in right cluster and
    ask agents [set group [who] of min-one-of centroids [opinion-distance]]

    ;; Computation of centroids possitions
    compute-centroids-positions
  ]

  ;; Killing centroids without connected agents
  ask centroids [
    let wom who
    if (not any? agents with [group = wom]) [
      set N_centroids N_centroids - 1
      die
    ]
  ]

  ;show count centroids
  ;; Catching the run-time error:
  ;; If there is just one component since all agents has same opinion, then the polarisation algorithm does produce error --
  ;; because of computing mean of empty list of distences: this list is empty since the only one existing centroid can't
  ;; compute distance to itself via double 'while' structuresince 'ai' and 'aj' lists are empty.
  ;; In this case it is obvious that polarisation is 0, so we set 'polarisation' and 'normalized_polarisation' to 0 directly via 'ifelse' structure
  ifelse (count centroids < 2) [
    ;show "Manual setting of polarisation globals to 0!"
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
    ;show whos
    ask agents [set distance_to_centroid [opinion-distance] of centroid group]  ;; Each agent computes her distance to her centroid and stores it as 'distance_to_centroid'.

    ;; Computing polarization -- distances of all centroids
    let ai but-last whos  ;; List of all 'i' -- 'who' initializing distances computation
    let aj but-first whos  ;; List of all 'j' -- 'who' of other end of distances computation
    foreach ai [i ->
      foreach aj [j ->
        ;show (word i ";  " j)
        ;; Each distance is weighted by fraction of both centroid groups
        ;; via formula '(N_centroids ^ 2) * (count agents with [group = i] / count agents) * (count agents with [group = j] / count agents)'
        let weight (N_centroids ^ 2) * (count agents with [group = i] / count agents) * (count agents with [group = j] / count agents)
        let cent-dist opinion-distance3 ([opinion-position] of centroid i) ([opinion-position] of centroid j)
        set distances lput (weight * cent-dist) distances
        set unweighted_distances lput cent-dist unweighted_distances
      ]
      set aj but-first aj
    ]

    ;; Computing polarization -- diversity in groups
    foreach sort [who] of centroids [wg ->
      let weight (count agents with [group = wg] / count agents)
      let cent-div (mean [distance_to_centroid] of agents with [group = wg])
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
    ;show diversity
    ;show distances
    ;show (word polarisation ";  " normalized_polarisation "; " unweighted_polarisation ";  " unweighted_normalized_polarisation)
  ]

  ;; Final coloring and killing of centroids
  if centroid_color? [ask agents [set color [color] of centroid group]]
  if killing_centroids? [ask centroids [die]]
  ;set N_centroids original_centroids_value
end


;; Sub-routine of polarization routine
to compute-centroids-positions
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
          set Opinion-position replace-item dim Opinion-position mean [item dim Opinion-position] of agents with [group = grp]
          set dim dim + 1
        ]
      ]
    ]
    set grp grp + 1
  ]
  ask centroids [
    getPlace
    ;show opinion-distance3 (Last-opinion) (Opinion-position)
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
  ;; Redundant conditions which should be avoided -- if the boundary is drawn as constant, then is completely same whether agents vaguely speak or openly listen,
  ;; seme case is for probability speaking of 100%, then it's same whether individual probability is drawn as constant or uniform, result is still same: all agents has probability 100%.
  if avoid-redundancies? and mode = "vaguely-speak" and boundary-drawn = "constant" [stop]
  if avoid-redundancies? and p-speaking-level = 1 and p-speaking-drawn = "uniform" [stop]
  ;; these two conditions cover 7/16 of all simulations, approx. the half! This code should stop them from running.

  ;; Check whether we set properly parameter 'updating' --
  ;; if we want update more dimensions than exists in simulation, then we set 'updating' to max of dimensions, i.e. 'opinions'
  if updating > opinions [set updating opinions]
end


;; Just envelope for updating agent at the begining of GO procedure
to preparing-myself
    set speak? speaking
    getColor
    getPlace

    ;; storing previous opinion position as 'Last-opinion'
    set Last-opinion Opinion-position

    ;; Firstly we have to determine dissatisfaction with the neighborhood
    set Satisfied? get-satisfaction
end


to set-group-identities
    ;; Preparation
  ask centroids [die]
  ;let original_centroids_value N_centroids

  ;; Detection of clusters via Louvain
  ask l-distances [die]  ;; Cleaning environment
  ask agents [create-l-distances-with other agents with [(sqrt(4 * opinions) - opinion-distance) / sqrt(4 * opinions) >= id_threshold] [set l-weight opinion-distance2 ([opinion-position] of end1)([opinion-position] of end2)]]
  ;show count l-distances
  ;ask l-distances with [l-weight < d_threshold] [die]
  nw:set-context agents l-distances ;with [l-weight >= d_threshold]
  let communities nw:louvain-communities
  ;show count l-distances
  ask l-distances [die]
  set N_centroids length communities

  ;; Computing clusters' mean 'opinion-position'
  let postions-clusters [] ;; List with all positions of all clusters
  foreach communities [c ->
    let one []  ;; List for one positio nof one cluster
    foreach range opinions [o ->
      set one lput precision (mean [item o opinion-position] of c) 3 one
    ]
    ;show one
    set postions-clusters lput one postions-clusters
    ;show postions-clusters
  ]

  ;; Preparation of centroids -- feedeing them with communities
  create-centroids N_centroids [
    set heading (who - min [who] of centroids)
    set Opinion-position item heading postions-clusters  ;; We set opinions, we try to do it smoothly...
    ;show Opinion-position
    set shape "circle"
    set size 1.5
    set color 5 + who * 10
    getPlace
  ]

  ;; Assignment of agents to groups
  ;let min_group min [who] of centroids
  ask agents [set group [who] of min-one-of centroids [opinion-distance]]

  ;; Computation of centroids possitions
  compute-centroids-positions

  ;let iter 0

;  ;; Iterating cycle -- looking for good match of centroids
;  while [sum [opinion-distance3 (Last-opinion) (Opinion-position)] of centroids > Centroids_change] [
;
;    ;set iter iter + 1
;    ;show (word "Iteration: " iter)
;
;    ;; turtles compute whether they are in right cluster and
;    ask agents [set group [who] of min-one-of centroids [opinion-distance]]
;
;    ;; Computation of centroids possitions
;    compute-centroids-positions
;  ]

  ;; Killing centroids without connected agents
  ask centroids [
    let wom who
    if (not any? agents with [group = wom]) [
      set N_centroids N_centroids - 1
      die
    ]
  ]

  ;; Final coloring and killing of centroids
  if centroid_color? [ask agents [set color [color] of centroid group]]
  if killing_centroids? [ask centroids [die]]
  ;set N_centroids original_centroids_value
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;   G O !   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; Main routine
to go
  ;; Just checking and avoiding runtime errors part of code
  avoiding-run-time-errors

  ;; Before a round erasing indicator of change
  set network-changes 0

  ;; Prepare group identities via Louvain
  set-group-identities

  ;; True part of GO procedure!
  ask agents [
    ;; speaking, coloring, updating and SATISFACTION!!!
    preparing-myself

    ;; Mechanism of own opinion or network change -- decision and rossolution
    ;; In case of dissatisfaction agents leave, otherwise updates opinion
    ifelse not satisfied? [
      change-of-network
      set network-changes network-changes + 1  ;; We advance the counter of network changes -- now, one just happened.
    ] [
      if model = "HK" [change-opinion-HK]
      ;; Note: Now here is only Hegselmann-Krause algorithm, but in the future we might easily employ other algorithms here!
    ]
  ]

  ;; The main algorithm might produce lonely agents, now we connect them to one other speaking agent
  connect-loners

  ;; Recoloring patches, computing how model settled down
  updating-patches-and-globals

  tick

  ;; Finishing condition:
  ;; 1) We reached state, where no turtle changes for RECORD-LENGTH steps, i.e. average of MAIN-RECORD (list of averages of turtles/agents RECORD) is 1 or
  ;; 2) We reached number of steps specified in MAX-TICKS
  if ((mean main-Record = 1 and network-changes <= 5) or ticks = max-ticks) [compute-polarisation-repeatedly]
  if ((mean main-Record = 1 and network-changes <= 5) or ticks = max-ticks) and record? [record-state-of-simulation]
  if (mean main-Record = 1 and network-changes <= 5) or ticks = max-ticks [stop]
  if (ticks / polarisation-each-n-steps) = floor (ticks / polarisation-each-n-steps) [compute-polarisation-repeatedly]
  if (ticks / record-each-n-steps) = floor(ticks / record-each-n-steps) and record? [record-state-of-simulation]
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
end


;; sub-routine computing whether the agent is satisfied with the neighborhood
to-report get-satisfaction
  ;; initialization of agent set 'supporters' containing agents whose opinion positions are key for agent's satisfaction
  let supporters nobody
  ;; 1) updating agent uses only visible link neighbors in small-world network
  let visibles comm-neighbors with [speak?]

  ;; 2) we have different modes for finding supporters:
  ;; 2.1) in mode "I got them!" agent looks inside her boundary (opinion +/- uncertainty),
  ;;      i.e. agent takes opinions not that much far from her opinion
  if mode = "openly-listen" [
    ;; we compute 'lim-dist' -- it is the numerical distance in given opinion space
    let lim-dist (Uncertainty * sqrt(opinions * 4))
    ;; we set as influentials agents with opinion not further than 'lim-dist'
    set supporters visibles with [opinion-distance <= lim-dist]
  ]

  ;; 2.1) in mode "They got me!" agent looks inside whose boundaries (opinion +/- uncertainty)
  ;;       she is, i.e. agents takes opinions spoken with such a big uncertainty that it matches her own opinion
  if mode = "vaguely-speak"  [
    ;; Note: Here is used the 'Uncetainty' value of called agent, agent who might be used for updating,
    ;;       not 'Uncertainty' of calling agent who updates her opinion.
    set supporters visibles with [opinion-distance <= (Uncertainty * sqrt(opinions * 4))]
  ]
  ;print ((1 - (count supporters / count visibles)) < Tolerance)
  ;; Now we can return the True/False value, whether the agent is satisfied and among visible neighbors are enough supporters
  report ifelse-value (count visibles > 0) [((1 - (count supporters / count visibles)) < Tolerance)][TRUE]  ;; In case no visibles are in the neighborhood, then agent is happy.
end


;; envelope controlling the way, how we change the network
to change-of-network
  if network-change = "link" [rewire-the-most-annoying-link]
  if network-change = "community" [leave-the-neighborhood-join-a-new-one]
  ;; Note: here might be other ways in the future, that's why the 'ifelse' structure is not used here
end


;; subroutine for leaving the neighborhood and joining a new one -- agent is decided to leave, we just process it here
to leave-the-neighborhood-join-a-new-one
  ;; Firstly, we have to count agents neighbors, to determine how many links agent has to create in the main part of the procedure
  let to-visibles my-comms with [(end1 != myself and [speak?] of end1) or (end2 != myself and [speak?] of end2)]
  let nei-size count to-visibles

  ;; Secondly, we cut off all the links
  ;show (word "I'm killing " nei-size " visible neighbors!")
  ask to-visibles [die]

  ;; Catching possible error with not enough visible agents for creating 'comms'
  let speaking-others other agents with [speak?]
  if (nei-size > count speaking-others) [set nei-size count speaking-others ;show "Not enough visibles!"
  ]

  ;; Thirdly, random VS intentional construction of new neighborhood.
  ifelse random-network-change? [
    ;; We set new neighborhood randomly or...
    create-comms-with n-of nei-size speaking-others
    ;show (word "I'm creating " nei-size " links with random visible neighbors!")
  ][
    ;; ...creates it out of the closest neighbors.
    create-comms-with min-n-of nei-size speaking-others [opinion-distance]
    ;show (word "I'm creating " nei-size " links with closest visible neighbors!")
  ]

  ;; P.S. Just hiding links for better speed of code -- when we change/cut a link, all links become visible and that slows down the simulation.
  ask comms [set hidden? TRUE]
end


;; TO-DO: agents should cut-off only neighbors that they previously heard speak,
;;        we probably should create their memory whom they heard speak and onlythose agents might cut-off.
;;
;; Note:  Now I implement it in modest variant: agent cuts off the most annoying presently speaking agent --
;;        there must be at least one, since agents are satisfied by the rule with the empty neighborhood and
;;        they update neighborhood only in case of dissatisfaction.
;; DONE!
;;

;; subroutine for changing one link
to rewire-the-most-annoying-link
  ;; Firstly, we cut the link with speaking agent with the most different opinion
  let visibles comm-neighbors with [speak?]
  ;show visibles
  ifelse random-network-change? [
    let a-visible one-of visibles
    ask one-of my-comms  with [other-end = a-visible] [;show self
      die]
    ;show (word "One random link to visible " a-visible " killed!")
  ][
    let annoyer max-one-of comm-neighbors with [speak?] [opinion-distance]
    ask one-of my-comms with [other-end = annoyer] [;show self
      die]
    ;show (word "Link to most annoying visible " annoyer " killed!")
  ]

  ;; Secondly, we choose for the agent a new speaking partner with the most close opinion
  let potentials other agents with [speak? and not comm-neighbor? myself]
  ;show potentials
  ifelse random-network-change? [
    create-comm-with one-of potentials ;[show self]
    ;show (word "Link to One random visible created!")
  ][
    let partner min-one-of potentials [opinion-distance]
    create-comm-with partner ;[show self]
    ;show (word "Link to the closest visible " partner " created!")
  ]

  ;; P.S. Just hiding links for better speed of code -- when we change/cut a link, all links become visible and that slows down the simulation.
  ask comms [set hidden? TRUE]
end


to connect-loners
  ;; We check whether each agent has at least one neighbor
  ask agents with [(count comm-neighbors) = 0] [
    ;; NOTE: Potential BUG! In case the agent without neis is the only speaking agent then 'potentials' = NOBODY and
    ;; it produces BUG during link creation.
    ;; That's why I catch it via 'if' structure -- if there is noone speaking, then the lone agent has to wait until the next round.
    ifelse (count other agents with [speak?] > 0) [
      let potentials other agents with [speak?]  ; We set 'potentials' to all other speaking agents and then...
      ;show "Creating new link!"
      create-comm-with ifelse-value (random-network-change?) [one-of potentials][min-one-of potentials [opinion-distance]] ;[show myself]  ;... it depends on scenario: we choose randomly or with the closest opinion
      ;print "Link just has been added!"
    ][;show "Not any speaking agents!"
    ]
  ]

  ;; P.S. Just hiding links for better speed of code -- when we change/cut a link, all links become visible and that slows down the simulation.
  ask comms [set hidden? TRUE]
end



;; sub-routine for updating opinion position of turtle according the Hegselmann-Krause (2002) model
to change-opinion-HK
  ;; initialization of agent set 'influentials' containing agents whose opinion positions uses updating agent
  let influentials nobody
  ;; 1) updating agent uses only visible link neighbors in small-world network
  let visibles other comm-neighbors with [speak?]

  ;; 2) we have different modes for finding influentials:
  ;; 2.1) in mode "I got them!" agent looks inside his boundary (opinion +/- uncertainty),
  ;;      i.e. agent takes opinions not that much far from her opinion
  if mode = "openly-listen" [
    ;; we compute 'lim-dist' -- it is the numerical distance in given opinion space
    let lim-dist (Uncertainty * sqrt(opinions * 4))
    ;; we set as influentials agents with opinion not further than 'lim-dist'
    set influentials visibles with [opinion-distance <= lim-dist]
  ]

  ;; 2.1) in mode "They got me!" agent looks inside whose boundaries (opinion +/- uncertainty)
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
  report precision weight 3
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
  report precision dist 8
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
182
75
N-agents
N-agents
10
1000
257.0
1
1
NIL
HORIZONTAL

SLIDER
10
385
102
418
n-neis
n-neis
1
500
46.0
1
1
NIL
HORIZONTAL

SLIDER
10
417
102
450
p-random
p-random
0
0.5
0.1
0.01
1
NIL
HORIZONTAL

SLIDER
10
75
182
108
opinions
opinions
1
50
2.0
1
1
NIL
HORIZONTAL

BUTTON
271
495
326
528
getPlace
ask turtles [getPlace]
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
140
102
185
model
model
"HK"
0

SLIDER
10
306
182
339
p-speaking-level
p-speaking-level
0
1
0.474
0.001
1
NIL
HORIZONTAL

SLIDER
10
228
182
261
boundary
boundary
0.01
1
0.2
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
-1.36807062E9
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
328
495
383
528
getColor
ask turtles [\n  set speak? TRUE\n  getColor\n]
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
381
495
436
528
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
436
495
491
528
sizes
show sort remove-duplicates [count turtles-here] of turtles 
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
271
462
326
495
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
325
462
380
495
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
261
148
306
boundary-drawn
boundary-drawn
"constant" "uniform"
1

PLOT
967
256
1167
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
380
462
468
495
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
467
462
529
495
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
529
462
589
495
Show links
ask links [set hidden? FALSE]
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
2.0
1
1
NIL
HORIZONTAL

BUTTON
588
462
659
495
avg. Opinion
show mean [mean opinion-position] of turtles
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

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
781
74
919
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
1001
227
Stability of turtles (average)
NIL
NIL
0.0
10.0
0.0
1.0
true
true
"" ""
PENS
"Turtles" 1.0 0 -16777216 true "" "plot mean [mean Record] of agents"
"Main record" 1.0 0 -2674135 true "" "plot mean main-Record"
"Normalized" 1.0 0 -13791810 true "" "plot normalized_polarisation"
"Polarisation" 1.0 0 -13345367 true "" "plot polarisation"

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
25.0
1
1
NIL
HORIZONTAL

CHOOSER
10
184
148
229
mode
mode
"openly-listen" "vaguely-speak"
0

PLOT
967
375
1167
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
-1368070620_257_0.1_46_2_1_0.2_uniform_0.474_uniform_openly-listen
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
0
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
2500.0
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
339
148
384
p-speaking-drawn
p-speaking-drawn
"constant" "uniform" "function"
1

PLOT
1166
375
1365
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
5000.0
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
107
182
140
tolerance-level
tolerance-level
0
1.1
0.5
0.01
1
NIL
HORIZONTAL

CHOOSER
101
140
193
185
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
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot network-changes / count turtles * 100"

SLIDER
9
495
130
528
conformity-level
conformity-level
0
1
0.5
0.01
1
NIL
HORIZONTAL

CHOOSER
130
495
249
540
conformity-drawn
conformity-drawn
"constant" "uniform"
1

CHOOSER
10
450
102
495
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
102
462
271
495
random-network-change?
random-network-change?
1
1
-1000

BUTTON
545
519
620
552
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
9.0
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
1166
256
1365
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
691
507
768
552
NIL
polarisation
17
1
11

MONITOR
967
507
1091
552
NIL
normalized_polarisation
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
0
1
-1000

SWITCH
1171
224
1319
257
killing_centroids?
killing_centroids?
1
1
-1000

SLIDER
1021
160
1171
193
id_threshold
id_threshold
0.01
1
0.9
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
1
10000
100.0
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
9
527
121
560
use_clusters?
use_clusters?
1
1
-1000

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
