;; Model for study of mutual interaction between public opinion and network structure.
;; RQ: How does identity and opinion openness influence polarisation of public opinion?
;;
;; This code is from project of Mike, Ashley, Ashwin and FranCesko,
;; we apply Hegselmann-Krausse model in more than 1D and look how agents adapt in >1D opinion space and whether they form groups,
;; then we include individual assignment of openness and identity sensitivity.
;;
;; !!!EXPERIMENTAL FEATURE: IDENTITY!!!
;; !!!SUPER EXPERIMENTAL: INDIVIDUAL IDENTITY VISION!!!
;;

;; Created:  2021-10-21 FranCesko
;; Edited:   2022-09-30 FranCesko
;; Encoding: windows-1250
;; NetLogo:  6.2.2
;;

;; IDEA: Choose, how many opinions agents update: sometimes 1, 2, 3, 4 ...
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
;;    b) separate file with agent's traits (P-speaking, HK Boundary etc.)
;;    c) as it was before, contextual variables of one whole simulation run are coded in the filenames
;; DONE!
;;
;; 6) Implement K-clusters algorithm for addressing just 2 clusters.
;; DONE!
;;
;; 7) compute centroids distances, store them in table of matrices for each id sensitivity level, store minimal WHO of centroid/group ID and number of groups for each level and
;;    use all of these to compute id sigmoid, instead of letting agents forming id groups as agentset selection in each identity check.
;; DONE!
;;
;; 8) debug code for the 'Use_Identity?' FALSE, now code works only with 'Use_Identity?' TRUE
;;
;;
;;

extensions [nw matrix table]

breed [centroids centroid]
undirected-link-breed [l-distances l-distance]

turtles-own [own-opinion own-previous-opinion own-boundary own-record own-conformity own-SPIRO own-WhichGroupHasEachSPIROSortedMeIn own-group-number own-distance-to-centroid
  own-opinion-sigmoid-xOffset own-opinion-sigmoid-steepness own-opinion-dice? own-identity-sigmoid-xOffset own-identity-sigmoid-steepness own-identity-dice?]
l-distances-own [l-weight]

globals [main-Record agents positions_clusters ESBG_polarisation SPIRO_set
         IDs-and-ns-of-id-groups distance-matrices]


;; Initialization and setup
to setup
  ;; We erase the world and clean patches
  ca
  ask patches [set pcolor patch-color]

  ;; To avoid some random artificialities we have to set random seed.
  ;; But for BehaviorSearch we have to avoid seed control, since BeahaviorSearch does it.
  if not avoid_seed_control? [  ;; Note: if switch 'avoid_seed_control?' is TRUE, then we avoid setting the seed.
    let seed ifelse-value (set-seed?) [RS][new-seed]
    random-seed seed
    if not set-seed? [set RS seed]
  ]

  ;; Then we migh initialize agents/turtles
  crt Number_Of_Agents [
    set own-opinion n-values Number_Of_Opinion_Dimensions [precision (1 - random-float 2) 3]  ;; We set opinions...
    set own-previous-opinion own-opinion  ;; ...set last opinion as present opinion...
    set own-record n-values record-length [0]  ;; ... we prepare indicator of turtle's stability, at all positions we set 0 as non-stability...
    set own-conformity get-conformity  ;; setting individual conformity level, and ...
    set own-boundary get-HK-boundary  ;;... setting value of HK boundary.
    set own-SPIRO get-own-SPIRO  ;; Individual sensitivity for group tightness/threshold.
    get-sigmoids ;;getting sigmoid parameters for opinion and identity influence probabilities
    getColor  ;; Coloring the agents according their opinion.
    getPlace  ;; Moving agents to the opinion space according their opinions.
  ]

  ;; We also have to set up some globals and create links for measuring opinion distances:
  set agents turtle-set turtles  ;; Note: If we just write 'set agents turtles', then variable 'agents' is a synonym for 'turtles', so it will contain in the future created centroids!

  ;; In case we don't use identity, we don't need to create links
  if Use_Identity? [
    ask agents [create-l-distances-with other agents]  ;; Creating full network for computing groups and polarisation
    update-l-distances-weights  ;; Setting distances links' weights
    ask l-distances [set hidden? true]  ;; Hiding links for saving comp. resources
  ]

  ;; Setting identity levels according identity scenario:
  ifelse Use_Identity? and Identity_Type = "individual" [
    ;; If we use 'individual' perceptions of identity groups, we have to set up levels of identity sensitivity:
    if Identity_Levels = 1 [set Maximum_SPIRO SPIRO_Mean]
    compute-identity-thresholds
  ][
    ;; If we use 'global' perception, then there is only one level and each agent has same value of 'own-SPIRO':
    ;; Firstly, we set 'SPIRO_set' as list of one constant value: 'SPIRO_Mean':
    set SPIRO_set (list SPIRO_Mean)

    ;; Secondly, we set 'own-SPIRO' of agents to the constant value:
    ask agents [
      set own-SPIRO SPIRO_Mean
    ]
  ]

  ;; Setting agents' identity groups
  ;; Everything is prepared for equal processing in all three cases of 'non-identity', 'global' and 'individual' perception of identity groups,
  ;; that's why we use only one procedure here for all three scenarios. But we handle them equally:
  ;; we process it for all identity levels -- in case of non-identity and global, there is just one.
  if Use_Identity? [set-group-identities]

  ;; Coloring patches according the number of agents/turtles on them.
  ask patches [set pcolor patch-color]

  ;; Setting the indicator of change for the whole simulation, again as non-stable.
  set Main-record n-values record-length [0]

  ;; Compute polarisation
  compute-polarisation-repeatedly

  reset-ticks
end


;; for computing thresholds and assigning levels to each agent as per drawing parameters
to compute-identity-thresholds
  ;; Computing values:
  ;; We know that SPIRO_Mean lower than 0.4 produces one identity group, that's why the most tolerant group/level
  ;; will have threshold 0.4. We also know that thresholds beyond 0.80 fracture the public to many groups and it makes
  ;; no sense use higher threshold than 0.8, so that's why the most intolerant group/level will have threshold 0.8.
  ;; If we will use just two identity levels, they will be 0.4 and 0.8. If we will use more levels, we will smoothly
  ;; distribute them between 0.4 and 0.8.
  set SPIRO_set n-values Identity_Levels [Maximum_SPIRO]
  foreach range (Identity_Levels - 1) [i -> set SPIRO_set replace-item i SPIRO_set precision (Minimum_SPIRO + (i * ((Maximum_SPIRO - Minimum_SPIRO) / (Identity_Levels - 1)))) 3]


  ;; We "ceiling" values 'own-SPIRO' of agents to the values of 'SPIRO_set', i.e.
  ;; we find the closest higher value of 'SPIRO_set' to the 'own-SPIRO' of agent and
  ;; change it for the closest value of 'SPIRO_set'.
  ask agents [
    ;; Firstly, we have to process the 'SPIRO_set': substract it from 'own-SPIRO' of agent.
    let diff map [v -> abs(v - own-SPIRO)] SPIRO_set

    ;; Secondly, we have to find the position of minimal 'diff':
    let pos position (min diff) diff

    ;; Thirdly, we have to check, whether the closest value of 'SPIRO_set' is higher,
    ;; if not, we have to point agent to the higher value, i.e. 'set pos pos  + 1",
    ;; and also check whether the 'pos' points on correct items on 'SPIRO_set', i.e.
    ;; whether we do/not jump out of range. In case we point 'outside the list',
    ;; we have to set 'pos' to maximal value.
    if not (pos = (length SPIRO_set - 1) or own-SPIRO < item pos SPIRO_set) [set pos pos + 1]

    ;; Finally, we set 'own-SPIRO' as the closest value of 'SPIRO_set':
    set own-SPIRO item pos SPIRO_set
  ]

  ;; Very lastly, we have to check whether all values of 'SPIRO_set' are represented in the agents,
  ;; it is possible, that for some low SD of random normal distributions of 'SPIRO_Mean' some values of
  ;; 'SPIRO_set' would not be represented by any agent.
  ;; Then it is obsolete to perform Louvain for non-represented values of 'SPIRO_set'. So...
  ;; Now we go through 'SPIRO_set' list value by value and only represented remain.
  set SPIRO_set remove-duplicates [own-SPIRO] of agents
end


; Setting identity groups via threshold levels to account for differences in sensitivity to group relationships
to set-group-identities
  ;; Firstly, we have to erase lists 'own-WhichGroupHasEachSPIROSortedMeIn' of agents,
  ;; since every step we have to fill it by new IDs of identity centroids:
  ask agents [set own-WhichGroupHasEachSPIROSortedMeIn table:make]
  ;; We erase/create these global tables, as well.
  set IDs-and-ns-of-id-groups table:make  ;; NOTE: the key is the level of group sensitivity, stored is the list: the first value in the list is number of groups, second/last is the lowest ID of group centroid.
  set distance-matrices table:make

  ;; Secondly, we go one effective ID threshold level after another, perform Louvain and k-means clusters for every level and store memberships.
  foreach SPIRO_set [ idtl ->
    ;; Cleaning environment
    ask centroids [die]

    ;; Detection of clusters via Louvain: Detection itself
    let selected-agents agents with [2 <= count my-l-distances with [l-weight >= idtl]]  ;; Note: We take into account only not loosely connected agents
    nw:set-context selected-agents l-distances with [l-weight >= idtl]  ;; For starting centroids we take into account only not loosely connected agents, but later we set groups for all.
    let communities nw:louvain-communities
    set N_centroids length communities

    ;; Computing clusters' mean 'own-opinion'
    set positions_clusters [] ;; List with all positions of all clusters
    foreach communities [c ->
      let one []  ;; List for one position of one cluster
      foreach range Number_Of_Opinion_Dimensions [o -> set one lput precision (mean [item o own-opinion] of c) 3 one]
      set positions_clusters lput one positions_clusters
    ]

    ;; Preparation of centroids -- feedeing them with communities
    create-centroids N_centroids [
      set heading (who - min [who] of centroids)
      set own-opinion item heading positions_clusters  ;; We set opinions, we try to do it smoothly...
      set shape "circle"
      set size 1.5
      set color 5 + (who - min [who] of centroids) * 10
      getPlace
    ]

    ;; Assignment of agents to groups
    ask agents [set own-group-number [who] of min-one-of centroids [opinion-distance]]  ;; Sic! Here we intentionally use all agents, including loosely connected.

    ;; Computation of centroids possitions
    compute-centroids-positions (agents)

    ;; Iterating cycle -- looking for good match of centroids
    while [sum [opinion-distance3 (own-previous-opinion) (own-opinion)] of centroids > Centroids_change] [

      ;; turtles compute whether they are in right cluster and
      ask agents [set own-group-number [who] of min-one-of centroids [opinion-distance]]

      ;; Computation of centroids possitions
      compute-centroids-positions (agents)
    ]

    ;; Storing ID of in-groups for the present identity level
    ask agents [table:put own-WhichGroupHasEachSPIROSortedMeIn idtl own-group-number]

    ;; Killing centroids without connected agents
    ask centroids [
      let wom who
      if (not any? agents with [own-group-number = wom]) [die]
    ]
    set N_centroids count centroids

    ;;;; Storing numbers and IDs of groups and centroid distances in tables
    ;; Firstly, number and lowest ID, it's the easiest
    table:put IDs-and-ns-of-id-groups idtl (list N_centroids min [who] of centroids)  ;; NOTE: the first value in the list is number of groups, second/last is the lowest ID of group centroid.

    ;; Secondly, we fill in respective item in 'distance-matrices'
      ;; We create distance matrix ...
      let min-id last table:get IDs-and-ns-of-id-groups idtl  ;; We have to know ID of the first centroid of respective level...
      let max-id (min-id - 1 + first table:get IDs-and-ns-of-id-groups idtl)  ;; ... so do we have to know the last centroid's ID.
      let m []  ;; future distance matrix -- initialized as empty row list
      let i min-id
      while [i <= max-id] [
        let j min-id
        let row []  ;; future row of distance matrix -- initialized as empty list
        while [j <= max-id] [
          set row lput ifelse-value (j = i) [0][precision (opinion-distance3 ([own-opinion] of centroid i) ([own-opinion] of centroid j)) 3] row
          set j j + 1
        ]
        set m lput row m
        set i i + 1
      ]

      ;; And then save it as the table entry...
      table:put distance-matrices idtl matrix:from-row-list m

    ;; Final coloring and killing of centroids
    if centroid_color? [ask agents [set color (5 + 10 * (own-group-number - min [who] of centroids))]]
    if killing_centroids? [ask centroids [die]]
  ]
  ;print IDs-and-ns-of-id-groups
  ;print distance-matrices
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;   G O !   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; Main routine
to go
  ;;;; Preparation part ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; All preparations of agents, globals, avoiding errors etc. in one sub-routine
  prepare-everything-for-the-step

  ;;;; Main part ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ask agents [
    if model = "HK" [change-opinion-HK]
    ;; Note: Now here is only Hegselmann-Krause algorithm, but in the future we might easily employ other algorithms here!
  ]

  ;;;; Final part ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Recoloring patches, agents, computing how model settled down
  updating-patches-and-globals

  tick

  ;; Finishing condition:
  ;; 1) We reached number of steps specified in MAX-TICKS
  recording-situation-and-computing-polarisation
  if ticks = max-ticks [stop]
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;   SUB/PROCEDURES   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to prepare-everything-for-the-step
  ;; Just checking and avoiding runtime errors part of code
  avoiding-run-time-errors

  ;; Update group identities via Louvain -- only if we use group identity.
  if Use_Identity? [
    update-l-distances-weights
    set-group-identities
  ]

  ;; Coloring and updating
  ask agents [preparing-myself]
end


to preparing-myself
  ;; Updating color and place
  getColor
  getPlace

  ;; storing previous opinion position as 'own-previous-opinion'
  set own-previous-opinion own-opinion
end


;; sub-routine for updating opinion position of turtle according the Hegselmann-Krause (2002) model
to change-opinion-HK
  ;; We define all other agents as NEIGHBS, i.e. potential INFLUENTIALS
  let neighbs other agents

  ;; In the first block of code we have to determine INFLUENTIALS -- the agents to whom the  updating agents listens to
  ;; Since rolling identity dice is computationally less demanding, we start with id dice:
  ;; If we use identity, we check distance of group centroids of NEIGHBS and filter them out.
  ;; In case we dont't use identity, them NEIGHBS are still all other agents.
  if Use_Identity? [
    ;; Identity check -- preparation of myself from global variables:
    let my-group-id table:get own-WhichGroupHasEachSPIROSortedMeIn own-SPIRO
    let distance-matrix table:get Distance-matrices own-SPIRO
    let group-info table:get IDs-and-ns-of-id-groups own-SPIRO  ;; NOTE: the first value in the list is number of groups, second/last is the lowest ID of group centroid.
    let my-i (my-group-id - last group-info)  ;; 'my-i' is distance matrix column we want to use, so we subtract group ID of myself from minimal group ID in respective level

    ;; Sigmoid code:
    ask other agents [
      ;; Firstly, each neighbor has to find 'her-j', i.e. distance matrix row for identity check:
      let her-group-id table:get own-WhichGroupHasEachSPIROSortedMeIn [own-SPIRO] of myself
      let her-j  (her-group-id - last group-info)   ;; 'her-j' is distance matrix row we want to use, so we subtract group ID of self from minimal group ID in respective level

      ;; Secondly, we roll the identity dice...
      let our-distance matrix:get distance-matrix my-i her-j  ;; After all the computations we finally get distance of centroids from distance matrix...

      ;;NOTE: The commented out part below can be uncommented if one wants to speed up simulation by hard-setting probabilities of
      ;;extreme distance values ( 0 or 1 ) to extreme values (p = 1 or p = 0). This would lead to differences in outcome especially
      ;;when the sigmoid steepnesses are low.

;      ifelse our-distance > 0 and our-distance < 1 [
        roll-identity-dice (our-distance)
      ;][
;        ifelse our-distance = 0[
;          set identity-dice? true
;        ][
;          set identity-dice? false
;        ]
;      ]

    ]

    ;; Now we set successful AGENTS as NEIGHBS:
    set neighbs other agents with [own-identity-dice?]
    if show_dice_rolls? [print count neighbs]
  ]

  ;; NOW we roll probabilistic dice based on opinion distance from influencer - if an agent is too far they are less likely to be heard this time.
  ;; first find out which agent is going to be heard this time.
  ask neighbs [
    let op-dist opinion-distance3 (own-opinion) ([own-opinion] of myself) ;; Note: We compute opinion distance ('op-dist') of NEIGHBS member and updating agent and ...
    roll-opinion-dice (precision op-dist 3)  ;; ... pass it rounded to 3 digits as the argument of ROLL-OPINION-DICE function.
  ]

  ;; now only follow them, the successful rollers and set them as INFLUENTIALS...
  let influentials neighbs with [own-opinion-dice?]
  if show_dice_rolls? [print count influentials]

  ;; 3) we also add the updating agent into 'influentials'
  set influentials (turtle-set self influentials)

  ;; we check whether there is someone else then calling/updating agent in the agent set 'influentials'
  if count influentials > 1 [

    ;; here we draw a list of dimensions which we will update:
    ;; by 'range opinions' we generate list of integers from '0' to 'opinions - 1',
    ;; by 'n-of updating' we randomly take 'updating' number of integers from this list
    ;; by 'shuffle' we randomize order of the resulting list
    let op-list shuffle n-of updating range Number_Of_Opinion_Dimensions

    ;; we initialize counter 'step'
    let step 0

    ;; we go through the while-loop 'updating' times:
    while [step < updating] [
      ;; we initialize/set index of updated opinion dimension according the items on the 'op-list',
      ;; note: since we use while-loop, we go through each item of the 'op-list', step by step, iteration by iteration.
      let i (item step op-list)

      ;; ad 1: averge position computation
      let val  precision (mean [item i own-opinion] of influentials) 3 ;; NOTE: H-K model really assumes that agent adopts immediatelly the 'consesual' position
      ;; ad 2: updating/weighting 'val' by 'Conformity' and own opinion
      let my item i own-opinion
      set val my + ((val - my) * own-conformity)
      ;; ad 3: assigning the value 'val'
      set own-opinion replace-item i own-opinion val

      ;; advancement of counter 'step'
      set step step + 1
    ]
  ]
end


to roll-opinion-dice [op-dist]
  ifelse opinion_sigmoid? [
    let probability-of-interaction compute-sigmoid (op-dist) ([own-opinion-sigmoid-xOffset] of myself) ([own-opinion-sigmoid-steepness] of myself)
    let dice (random-float 1)
    set own-opinion-dice? dice < probability-of-interaction
    if show_dice_rolls? and probability-of-interaction > 0 and probability-of-interaction < 1 [
      print (word "Opinion distance: " op-dist ", Boundary: " [own-boundary] of myself ", Sigmoidial probability: " probability-of-interaction ", Rolled dice: " dice "; Result: " own-opinion-dice?)]
  ][set own-opinion-dice? op-dist <= [own-opinion-sigmoid-xOffset] of myself]  ;; If we don't use opinion sigmoid, then we set 'opinion-dice?' according sharp difference
end



to roll-identity-dice [our-distance]
  ifelse identity_sigmoid? [
    let probability-of-interaction compute-sigmoid (our-distance) ([own-identity-sigmoid-xOffset] of myself) ([own-identity-sigmoid-steepness] of myself)
    let dice (random-float 1)
    set own-identity-dice? dice < probability-of-interaction
    if show_dice_rolls? and probability-of-interaction > 0 and probability-of-interaction < 1 [
      print (word "Group distance: " our-distance ", Identity distance threshold: " [own-identity-sigmoid-xOffset] of myself ", Sigmoidial probability: " probability-of-interaction " Rolled dice: " dice "; Result: " own-identity-dice?)]
  ][set own-identity-dice? our-distance <= [own-identity-sigmoid-xOffset] of myself]  ;; If we don't use identity sigmoid, then we set 'identity-dice?' according sharp difference
end


;;subroutine for computing inverted (y decreases as x increases) sigmoid of input
to-report compute-sigmoid [x xOffset steepness]
  ;; Slider 'maxSteepness' determines maximum possible steepness of the sigmoids.
  ;; Sigmoid outputs from 1 to 0 as x increases from 0 to 1
  report precision (1 / (1 + exp (steepness * maxSteepness * (x - xOffset)))) 3
end


;; sub-routine for computing opinion distance of two comparing agents
to-report opinion-distance
  ;; we store in temporary variable the opinion of the called and compared agent
  let my own-opinion

  ;; we store in temporary variable the opinion of the calling and comparing agent
  let her [own-opinion] of myself

  ;; we initialize counter of step of comparison -- we will compare as many times as we have dimensions
  let step 0

  ;; we initialize container where we will store squared distance in each dimension
  let dist 0

  ;; while loop going through each dimension, computiong distance in each dimension, squarring it and adding in the container
  while [step < Number_Of_Opinion_Dimensions] [
    ;; computiong distance in each dimension, squarring it and adding in the container
    set dist dist + ((item step my - item step her) ^ 2)

    ;; advancing 'step' counter by 1
    set step step + 1
  ]

  ;; computing square-root of the container 'dist' -- computing Euclidean distance -- and setting it as 'dist'
  set dist sqrt dist

  ;; Computing normalization constant according switch 'normalize_sigmoid_distances?':
  let normalization ifelse-value (Normalize_Distances?) [1 / sqrt(4 * Number_Of_Opinion_Dimensions)][1]

  ;; reporting Euclidean distance
  report dist * normalization
end


;; sub-routine for computing opinion distance of two comparing opinion positions  -- relative distance weighted as 1 for minimal distance and 0 for the maximal one
to-report opinion-distance2 [my her]
  ;; we initialize counter of step of comparison -- we will compare as many times as we have dimensions
  let step 0

  ;; we initialize container where we will store squared distance in each dimension
  let dist 0

  ;; while loop going through each dimension, computiong distance in each dimension, squarring it and adding in the container
  while [step < Number_Of_Opinion_Dimensions] [
    ;; computiong distance in each dimension, squarring it and adding in the container
    set dist dist + (item step my - item step her) ^ 2

    ;; advancing 'step' counter by 1
    set step step + 1
  ]

  ;; computing square-root of the container 'dist' -- computing Euclidean distance -- and setting it as 'dist'
  set dist sqrt dist

  ;; Turning 'dist' into 'weight'
  let weight (sqrt(4 * Number_Of_Opinion_Dimensions) - dist) / sqrt(4 * Number_Of_Opinion_Dimensions)

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
  while [step < Number_Of_Opinion_Dimensions] [
    ;; computiong distance in each dimension, squarring it and adding in the container
    set dist dist + (item step my - item step her) ^ 2

    ;; advancing 'step' counter by 1
    set step step + 1
  ]

  ;; computing square-root of the container 'dist' -- computing Euclidean distance -- and setting it as 'dist'
  set dist sqrt dist

  ;; Computing normalization constant according switch 'normalize_sigmoid_distances?':
  let normalization ifelse-value (Normalize_Distances?) [1 / sqrt(4 * Number_Of_Opinion_Dimensions)][1]

  ;; reporting weight of distance
  report precision (dist * normalization) 10
end


to recording-situation-and-computing-polarisation
  ;; Recording condition:
  ;; 1) We reached end, e.g. number of steps specified in MAX-TICKS
  ;; 2) Recording and computing polarisation on the fly, e.g. we reached POLARISATION-EACH-N-STEPS
  if (ticks = max-ticks) or (ticks / polarisation-each-n-steps) = floor (ticks / polarisation-each-n-steps) [compute-polarisation-repeatedly]
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
    set own-record fput ifelse-value (own-previous-opinion = own-opinion) [1][0] but-last own-record

  ]
  ;; Then we might update it for the whole:
  set main-Record fput precision (mean [mean own-record] of agents) 3 but-last main-Record

  ;; Coloring agents according identity group
  if centroid_color? [ask agents [set color (5 + 10 * (own-group-number - min [own-group-number] of agents))]]
end


;; Procedure reporting ESBG/Ashwin's polarisation
to-report Ash-polarisation
  ;; Preparation
  create-centroids 2 [set shape "square" set own-opinion n-values Number_Of_Opinion_Dimensions [precision (1 - random-float 2) 8]]
  let cent1 max [who] of centroids  ;; Storing 'who' of two new centroids
  let cent0 cent1 - 1
  ;; Random assignment of agents to the groups -- we create random list of 'cent0s' and 'cent1s' of same length as agents number and then assign them based on agent's WHO
  let membership shuffle (sentence n-values round (Number_Of_Agents / 2) [cent0] n-values (ifelse-value (0 = Number_Of_Agents mod 2) [Number_Of_Agents / 2][(Number_Of_Agents - 1) / 2]) [cent1])
  ask agents [set own-group-number item who membership]
  updating-centroids-own-opinion (cent0) (cent1)  ;; Initial update

  ;; Iterating until centroids are stable
  while [Centroids_change < sum [opinion-distance3 (own-opinion) (own-previous-opinion)] of centroids with [who >= cent0]][
    update-agents-opinion-group (cent0) (cent1)
    updating-centroids-own-opinion (cent0) (cent1)
  ]

  ;; Computing polarisation -- cutting-out agents too distant from centroids
  ask agents [set own-distance-to-centroid [opinion-distance] of centroid own-group-number]
  let a0 agents with [own-group-number = cent0]
  set a0 min-n-of (count a0 - ESBG_furthest_out) a0 [opinion-distance3 ([own-opinion] of self) ([own-opinion] of centroid cent0)]
  let a1 agents with [own-group-number = cent1]
  set a1 min-n-of (count a1 - ESBG_furthest_out) a1 [opinion-distance3 ([own-opinion] of self) ([own-opinion] of centroid cent1)]

  ;; Updating centroids and agents opinion position (without furthest agents)
  ask centroid cent0 [foreach range Number_Of_Opinion_Dimensions [o -> set own-opinion replace-item o own-opinion precision (mean [item o own-opinion] of a0) 8] getPlace]
  ask a0 [set own-distance-to-centroid [opinion-distance] of centroid cent0]
  ask centroid cent1 [foreach range Number_Of_Opinion_Dimensions [o -> set own-opinion replace-item o own-opinion precision (mean [item o own-opinion] of a1) 8] getPlace]
  ask a1 [set own-distance-to-centroid [opinion-distance] of centroid cent1]

  ;; Preparing final distances and diversity
  let normalization ifelse-value (Normalize_Distances?) [1][1 / sqrt(4 * Number_Of_Opinion_Dimensions)]
  ;; NOTE: If the switch 'Normalize_Distances?' is true, then functions 'opinion-distance' and 'opinion-distance3' compute normalized distances,
  ;; then we must avoid normalization here, but if the switch is false, then the function computes plain distance and we must normalize here --
  ;; so we use same concept and almost same code as in the functions 'opinion-distance' and 'opinion-distance3' here, but reversed:
  ;; the true switch means no normalization here, the false switch means normalization here, so in the end we normalize values in ESBG polarization exactly once.
  let cent-dist normalization * opinion-distance3 ([own-opinion] of centroid cent0) ([own-opinion] of centroid cent1)
  let div0 normalization * (mean [own-distance-to-centroid] of a0)
  let div1 normalization * (mean [own-distance-to-centroid] of a1)

  ;; Cleaning and reporting
  ask centroids with [who >= cent0] [die]
  report (cent-dist / (1 + div0 + div1))
end


to update-agents-opinion-group [cent0 cent1]
  ;; Checking the assignment -- is the assigned centroid the nearest? If not, reassign!
  ask agents [set own-group-number own-group-number - ([who] of min-one-of centroids with [who >= cent0] [opinion-distance])]  ; set color 15 + group * 10]
  let wrongly-at-grp0 turtle-set agents with [own-group-number = -1]  ;; they are in 0, but should be in 1: 0 - 1 = -1
  let wrongly-at-grp1 turtle-set agents with [own-group-number = 1]  ;; they are in 1, but should be in 0: 1 - 0 = 1
  ifelse count wrongly-at-grp0 = count wrongly-at-grp1 [
    ask agents [set own-group-number [who] of min-one-of centroids with [who >= cent0] [opinion-distance]]
  ][
    let peleton agents with [own-group-number = 0]
    ifelse count wrongly-at-grp0 < count wrongly-at-grp1 [
      set peleton (turtle-set peleton wrongly-at-grp0 max-n-of (count wrongly-at-grp0) wrongly-at-grp1 [opinion-distance3 ([own-opinion] of self) ([own-opinion] of centroid cent0)]) ;; all agents assigned correctly + smaller group of wrong + from bigger group 'n of size of smaller group'
      let stayed agents with [not member? self peleton]
      ask peleton [set own-group-number [who] of min-one-of centroids with [who >= cent0] [opinion-distance] ;set color 15 + group * 10
      ]
      ask stayed [set own-group-number cent1 ;set color 15 + group * 10
      ]
     ][
      set peleton (turtle-set peleton wrongly-at-grp1 max-n-of (count wrongly-at-grp1) wrongly-at-grp0 [opinion-distance3 ([own-opinion] of self) ([own-opinion] of centroid cent1)]) ;; all agents assigned correctly + smaller group of wrong + from bigger group 'n of size of smaller group'
      let stayed agents with [not member? self peleton]
      ask peleton [set own-group-number [who] of min-one-of centroids with [who >= cent0] [opinion-distance] ;set color 15 + group * 10
      ]
      ask stayed [set own-group-number cent0 ;set color 15 + group * 10
      ]
    ]
  ]
end


to updating-centroids-own-opinion [cent0 cent1]
  ;; Storing opinion as own-previous-opinion
  ask centroids with [who >= cent0] [set own-previous-opinion own-opinion]

  ;; Computing groups mean 'own-opinion'
  set positions_clusters [] ;; List with all positions of both 2 groups
  foreach range 2 [c ->
    let one-position []  ;; List for one position of one cluster
    foreach range Number_Of_Opinion_Dimensions [o -> set one-position lput precision (mean [item o own-opinion] of agents with [own-group-number = cent0 + c]) 8 one-position]
    set positions_clusters lput one-position positions_clusters
  ]

  ;; Setting opinions of centroids
  ask centroid cent0 [set own-opinion item 0 positions_clusters getPlace]
  ask centroid cent1 [set own-opinion item 1 positions_clusters getPlace]
end


;; Sub-routine for updating Distances links' weights,
;; according opinion distance of both their ends
to update-l-distances-weights
  ;; We use function 'opinion-distance2', which needs two opinion positions as input and
  ;; receives their distance as output, but this distance is converted to weight:
  ;; weight = 1 means that both positions are same, weight = 0 means that their distance is maximal,
  ;; i.e. both positions are in oposit corners of respective N-dimensional space.
  ask l-distances [set l-weight opinion-distance2 ([own-opinion] of end1) ([own-opinion] of end2)]
  ask l-distances [set hidden? TRUE]
end


;; We compute polarisation several times and then set it for the average
to compute-polarisation-repeatedly
  ;; Initialization of temporal variables
  let r 0
  let ap []
  update-l-distances-weights

  ;; Repeating cycle
  while [r < polar_repeats] [
    set ap lput Ash-polarisation ap
    set r r + 1
  ]

  ;; Setting variables back
  set ESBG_polarisation precision (mean ap) 3
end


;; Sub-routine of polarization routine
to compute-centroids-positions [sel-agents]
  ;; Preparation
  ask centroids [set own-previous-opinion own-opinion]

  ;; Computation of centroids possitions
  let grp min [who] of centroids
  while [grp <= (max [who] of centroids)] [
    ask centroid grp [
      ifelse (not any? agents with [own-group-number = grp]) [
        set own-opinion own-previous-opinion
      ][
        let dim 0
        while [dim < Number_Of_Opinion_Dimensions] [
          set own-opinion replace-item dim own-opinion mean [item dim own-opinion] of sel-agents with [own-group-number = grp]
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


;; Sub-routine for assigning value of conformity
to-report get-conformity
  ;; We have to initialize empty temporary variable
  let cValue 0

  ;; Then we draw the value according the chosen method
  if Conformity_Distribution = "constant" [set cValue Conformity_Mean + random-float 0]  ;; NOTE! 'random-float 0' is here for consuming one pseudorandom number to cunsume same number of pseudorandom numbers as "uniform
  if Conformity_Distribution = "uniform" [set cValue ifelse-value (Conformity_Mean < 0.5)
                                                           [precision (random-float (2 * Conformity_Mean)) 3]
                                                           [precision (1 - (random-float (2 * (1 - Conformity_Mean)))) 3]]
  if Conformity_Distribution = "normal" [ set cValue precision (random-normal Conformity_Mean Conformity_STD) 3
                                   while [cValue > 1 or cValue <= 0] [ set cValue precision (random-normal Conformity_Mean Conformity_STD) 3]
  ]
  report cValue
end


;; Sub-routine for assigning value of
to-report get-own-SPIRO
  ;; We have to initialize empty temporary variable
  let gtValue 0

  ;; Then we draw the value according the chosen method
  if SPIRO_Distribution = "constant" [set gtValue SPIRO_Mean + random-float 0]  ;; NOTE! 'random-float 0' is here for consuming one pseudorandom number to cunsume same number of pseudorandom numbers as "uniform
  if SPIRO_Distribution = "uniform" [set gtValue ifelse-value (SPIRO_Mean < 0.5)
                                                           [precision (random-float (2 * SPIRO_Mean)) 3]
                                                           [precision (1 - (random-float (2 * (1 - SPIRO_Mean)))) 3]]
  if SPIRO_Distribution = "normal" [ set gtValue precision (random-normal SPIRO_Mean SPIRO_STD) 3
                                  while [gtValue > 1 or gtValue <= 0] [ set gtValue precision (random-normal SPIRO_Mean SPIRO_STD) 3]
  ]
  report gtValue
end


;; sub-routine for assigning value of boundary to the agent
to-report get-HK-boundary
  ;; We have to initialize empty temporary variable
  let uValue 0

  ;; Then we draw the value according the chosen method
  if Boundary_Distribution = "constant" [set uValue Boundary_Mean + random-float 0]  ;; NOTE! 'random-float 0' is here for consuming one pseudorandom number to cunsume same number of pseudorandom numbers as "uniform"
  if Boundary_Distribution = "uniform" [set uValue precision (random-float (2 * Boundary_Mean)) 3]
  if Boundary_Distribution = "normal" [set uValue precision (random-normal Boundary_Mean Boundary_STD) 3
                                while [uValue > 1 or uValue <= 0] [ set uValue precision (random-normal  Boundary_Mean Boundary_STD) 3]
  ]
  ;; reporting value back for assigning
  report uValue
end


;; Sub procedure for generating needed parameters for sigmoids
to get-sigmoids
  ;; Opinion ones
  set own-opinion-sigmoid-xOffset own-boundary
  set own-opinion-sigmoid-steepness precision (random-normal Mean_Opinion_Sigmoid_Steepness STD_Opinion_Sigmoid_Steepness) 3
  while [own-opinion-sigmoid-steepness > 1 or own-opinion-sigmoid-steepness < 0] [
    set own-opinion-sigmoid-steepness precision (random-normal Mean_Opinion_Sigmoid_Steepness STD_Opinion_Sigmoid_Steepness) 3
  ]

  ;; Identity ones
  set own-identity-sigmoid-xOffset precision (random-normal Mean_Identity_Sigmoid_xOffset STD_Identity_Sigmoid_xOffset) 3
  while [own-identity-sigmoid-xOffset > 1 or own-identity-sigmoid-xOffset < 0] [
    set own-identity-sigmoid-xOffset precision (random-normal Mean_Identity_Sigmoid_xOffset STD_Identity_Sigmoid_xOffset) 3
  ]
  set own-identity-sigmoid-steepness precision (random-normal Mean_Identity_Sigmoid_Steepness STD_Identity_Sigmoid_Steepness) 3
  while [own-identity-sigmoid-steepness > 1 or own-identity-sigmoid-steepness < 0] [
    set own-identity-sigmoid-steepness precision (random-normal Mean_Identity_Sigmoid_Steepness STD_Identity_Sigmoid_Steepness) 3
  ]
end


;; sub-routine for graphical representation -- it takes two opinion dimension and gives the agent on XY coordinates accordingly
to getPlace
  ;; check whether our cosen dimension is not bigger than maximum of dimensions in the simulation
  if X-opinion > Number_Of_Opinion_Dimensions [set X-opinion 1]
  if Y-opinion > Number_Of_Opinion_Dimensions [set Y-opinion 1]

  ;; then we rotate the agent towards the future place
  facexy ((item (X-opinion - 1) own-opinion) * max-pxcor) ((item (Y-opinion - 1) own-opinion) * max-pycor)

  ;; lastly we move agent on the place given by opinion dimensions chosen for X and Y coordinates
  set xcor (item (X-opinion - 1) own-opinion) * max-pxcor
  set ycor (item (Y-opinion - 1) own-opinion) * max-pycor
end


;; sub routine for coloring agents according their average opinion across all dimensions --
;; useful for distinguishing agents with same displayed coordinates, but differing in other opinion dimensions,
;; then we see at one place agents with different colors.
to getColor
  ;; speaking agents are colored from very dark red (average -1) through red (average 0) to very light red (average +1)
  set color 15 + 4 * mean(own-opinion)
  set size (max-pxcor / 10)
end


;; sub-routine for visual purposes -- colors empty patches white, patches with some agents light green, with many agents dark green, with all agents black
to-report patch-color
  report 59.9 - (9.8 * (ln(1 + count turtles-here) / ln(Number_Of_Agents)))
end


;; Sub routine just for catching run-time errors
to avoiding-run-time-errors
  ;; Check whether we set properly parameter 'updating' --
  ;; if we want update more dimensions than exists in simulation, then we set 'updating' to max of dimensions, i.e. 'opinions'
  if updating > Number_Of_Opinion_Dimensions [set updating Number_Of_Opinion_Dimensions]
end
@#$#@#$#@
GRAPHICS-WINDOW
343
10
791
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
46
146
79
Number_Of_Agents
Number_Of_Agents
10
1000
129.0
1
1
NIL
HORIZONTAL

SLIDER
11
78
234
111
Number_Of_Opinion_Dimensions
Number_Of_Opinion_Dimensions
1
50
2.0
1
1
NIL
HORIZONTAL

CHOOSER
251
507
343
552
model
model
"HK"
0

SLIDER
1
276
130
309
Boundary_Mean
Boundary_Mean
0.0
1
0.2
0.001
1
NIL
HORIZONTAL

INPUTBOX
798
10
905
70
RS
1.648985926E9
1
0
Number

SWITCH
905
10
998
43
set-seed?
set-seed?
0
1
-1000

BUTTON
559
491
614
524
inspect
inspect turtle 0\nask turtle 0 [print own-opinion]
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
614
491
669
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
503
491
558
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
503
458
558
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
188
157
233
Boundary_Distribution
Boundary_Distribution
"constant" "uniform" "normal"
2

PLOT
1122
226
1282
346
 'Boundary' Distribution
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
"default" 0.05 1 -16777216 true "" "histogram [own-boundary] of agents"

BUTTON
558
458
671
491
avg. Boundary
show mean [own-boundary] of turtles
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
645
458
707
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
707
458
767
491
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
1055
10
1147
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
1055
42
1147
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

PLOT
798
195
1092
403
Development of opinions
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
"Op01" 1.0 0 -16777216 true "" "plot mean [item 0 own-opinion] of agents"
"Op02" 1.0 0 -7500403 true "" "if Number_Of_Opinion_Dimensions >= 2 [plot mean [item 1 own-opinion] of agents]"
"Op03" 1.0 0 -2674135 true "" "if Number_Of_Opinion_Dimensions >= 3 [plot mean [item 2 own-opinion] of agents]"
"Op04" 1.0 0 -955883 true "" "if Number_Of_Opinion_Dimensions >= 4 [plot mean [item 3 own-opinion] of agents]"
"Op05" 1.0 0 -6459832 true "" "if Number_Of_Opinion_Dimensions >= 5 [plot mean [item 4 own-opinion] of turtles]"
"Op06" 1.0 0 -1184463 true "" "if Number_Of_Opinion_Dimensions >= 6 [plot mean [item 5 own-opinion] of turtles]"
"Op07" 1.0 0 -10899396 true "" "if Number_Of_Opinion_Dimensions >= 7 [plot mean [item 6 own-opinion] of turtles]"
"Op08" 1.0 0 -13840069 true "" "if Number_Of_Opinion_Dimensions >= 8 [plot mean [item 7 own-opinion] of turtles]"
"Op09" 1.0 0 -14835848 true "" "if Number_Of_Opinion_Dimensions >= 9 [plot mean [item 8 own-opinion] of turtles]"
"Op10" 1.0 0 -11221820 true "" "if Number_Of_Opinion_Dimensions >= 10 [plot mean [item 9 own-opinion] of turtles]"
"Op11" 1.0 0 -13791810 true "" "if Number_Of_Opinion_Dimensions >= 11 [plot mean [item 10 own-opinion] of turtles]"

SLIDER
1055
74
1147
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
798
76
1056
196
Stability of turtles (average)
NIL
NIL
0.0
10.0
0.0
0.2
true
true
"" ""
PENS
"Turtles" 1.0 0 -16777216 true "" "plot mean [mean own-record] of agents"
"Main record" 1.0 0 -2674135 true "" "plot mean main-Record"
"ESBG" 1.0 0 -11221820 true "" "plot ESBG_polarisation"

MONITOR
1036
406
1108
451
avg. Record
mean [mean own-record] of agents
3
1
11

SLIDER
1147
10
1258
43
record-length
record-length
10
100
15.0
1
1
NIL
HORIZONTAL

SLIDER
1147
42
1258
75
max-ticks
max-ticks
100
10000
730.0
5
1
NIL
HORIZONTAL

SLIDER
1
355
142
388
Conformity_Mean
Conformity_Mean
0
1
0.563
0.001
1
NIL
HORIZONTAL

CHOOSER
1
310
156
355
Conformity_Distribution
Conformity_Distribution
"constant" "uniform" "normal"
2

BUTTON
669
491
744
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
429
461
502
521
N_centroids
1.0
1
0
Number

SLIDER
1147
74
1297
107
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
1122
106
1282
226
'Conformity' Distribution
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
"default" 0.05 1 -16777216 true "" "histogram [own-conformity] of agents"

MONITOR
933
406
1034
451
NIL
ESBG_polarisation
17
1
11

SWITCH
926
452
1047
485
centroid_color?
centroid_color?
0
1
-1000

SWITCH
798
439
927
472
killing_centroids?
killing_centroids?
0
1
-1000

SLIDER
10
233
128
266
SPIRO_Mean
SPIRO_Mean
0.01
1
0.44
0.001
1
NIL
HORIZONTAL

SLIDER
1259
10
1429
43
polarisation-each-n-steps
polarisation-each-n-steps
0
10000
400.0
10
1
NIL
HORIZONTAL

SLIDER
1258
41
1363
74
polar_repeats
polar_repeats
1
100
50.0
1
1
NIL
HORIZONTAL

SWITCH
1
587
116
620
use_identity?
use_identity?
0
1
-1000

SLIDER
823
471
927
504
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
798
405
922
438
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
157
188
264
233
SPIRO_Distribution
SPIRO_Distribution
"constant" "uniform" "normal"
2

PLOT
1282
106
1442
226
'id_threshold' Distribution
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
"default" 0.01 1 -16777216 true "" "histogram [own-SPIRO] of agents"

CHOOSER
11
143
103
188
Identity_Type
Identity_Type
"global" "individual"
1

SLIDER
101
143
207
176
Identity_Levels
Identity_Levels
1
10
7.0
1
1
NIL
HORIZONTAL

SLIDER
128
233
239
266
SPIRO_STD
SPIRO_STD
0
1
0.075
0.001
1
NIL
HORIZONTAL

SLIDER
173
354
313
387
Conformity_STD
Conformity_STD
0
1
0.222
0.001
1
NIL
HORIZONTAL

SLIDER
130
276
252
309
Boundary_STD
Boundary_STD
0
1
0.05
0.001
1
NIL
HORIZONTAL

SLIDER
1
387
247
420
Mean_Opinion_Sigmoid_Steepness
Mean_Opinion_Sigmoid_Steepness
0
1
0.5
0.01
1
NIL
HORIZONTAL

SLIDER
251
551
402
584
maxSteepness
maxSteepness
0
700
700.0
10
1
NIL
HORIZONTAL

SLIDER
2
420
246
453
STD_Opinion_Sigmoid_Steepness
STD_Opinion_Sigmoid_Steepness
0
1
0.1
0.001
1
NIL
HORIZONTAL

SLIDER
2
453
257
486
Mean_Identity_Sigmoid_Steepness
Mean_Identity_Sigmoid_Steepness
0
1
0.7
0.001
1
NIL
HORIZONTAL

SLIDER
2
486
249
519
STD_Identity_Sigmoid_Steepness
STD_Identity_Sigmoid_Steepness
0
1
0.125
0.001
1
NIL
HORIZONTAL

SLIDER
210
112
342
145
Minimum_SPIRO
Minimum_SPIRO
0
0.5
0.29
0.001
1
NIL
HORIZONTAL

SLIDER
210
145
342
178
Maximum_SPIRO
Maximum_SPIRO
0.55
1
0.89
0.001
1
NIL
HORIZONTAL

SWITCH
11
111
185
144
Normalize_Distances?
Normalize_Distances?
0
1
-1000

SWITCH
799
507
926
540
show_dice_rolls?
show_dice_rolls?
1
1
-1000

SLIDER
3
519
243
552
Mean_Identity_Sigmoid_xOffset
Mean_Identity_Sigmoid_xOffset
0
1
0.15
0.001
1
NIL
HORIZONTAL

SLIDER
3
552
235
585
STD_Identity_Sigmoid_xOffset
STD_Identity_Sigmoid_xOffset
0
1
0.1
0.001
1
NIL
HORIZONTAL

PLOT
1121
345
1282
465
'opinion-steepness' Distribution
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
"default" 0.05 1 -16777216 true "" "histogram [own-opinion-sigmoid-steepness] of agents "

PLOT
1282
345
1442
465
'identity-steepness' Distribution
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
"default" 0.05 1 -16777216 true "" "histogram [own-identity-sigmoid-steepness] of agents "

PLOT
1282
225
1442
345
'identity-xOffset' Distribution
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
"default" 0.05 1 -16777216 true "" "histogram [own-identity-sigmoid-xOffset] of agents "

SWITCH
199
387
329
420
opinion_sigmoid?
opinion_sigmoid?
0
1
-1000

SWITCH
204
453
333
486
identity_sigmoid?
identity_sigmoid?
0
1
-1000

SWITCH
905
43
1055
76
avoid_seed_control?
avoid_seed_control?
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
  <experiment name="RS001-200_complementOfHillClimbingSearchInClassicHK" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="366"/>
    <metric>ESBG_polarisation</metric>
    <steppedValueSet variable="RS" first="1" step="1" last="200"/>
    <enumeratedValueSet variable="set-seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="avoid_seed_control?">
      <value value="false"/>
    </enumeratedValueSet>
    <steppedValueSet variable="boundary" first="0.1" step="0.01" last="0.5"/>
    <steppedValueSet variable="conformity-level" first="0.1" step="0.1" last="1"/>
    <enumeratedValueSet variable="boundary-drawn">
      <value value="&quot;constant&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conformity-drawn">
      <value value="&quot;constant&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N-agents">
      <value value="129"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="opinions">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="updating">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="polar_repeats">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="id_threshold-std">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Y-opinion">
      <value value="1"/>
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
    <enumeratedValueSet variable="mean-identity-sigmoid-steepness">
      <value value="0.85"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-identity-sigmoid-steepness">
      <value value="0.35"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="opinion_sigmoid?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="polarisation-each-n-steps">
      <value value="400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show_dice_rolls?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="d_threshold">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ESBG_furthest_out">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="identity_sigmoid?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxSteepness">
      <value value="700"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="identity_levels">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-length">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="id_threshold">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="centroid_color?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="use_identity?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-opinion-sigmoid-steepness">
      <value value="0.45"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N_centroids">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conformity-std">
      <value value="0.15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-identity-sigmoid-xOffset">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="killing_centroids?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary-std">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-identity-sigmoid-xOffset">
      <value value="0.35"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_id_level">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_id_level">
      <value value="0.85"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="normalize_distances?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std-opinion-sigmoid-steepness">
      <value value="0.15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold_drawn">
      <value value="&quot;normal&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="identity_type">
      <value value="&quot;individual&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X-opinion">
      <value value="1"/>
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
