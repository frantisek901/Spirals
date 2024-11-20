extensions [ pathdir ] ;; owl ]

globals
[
  ERROR?
  test-level
  max-memos-received
  max-memory-size-used
  max-memory-size-current
  min-memory-size-current
  multi-run?
  directory-name
  logfile-name
  agentfile-name
  sourcefile-name
  run-identifier
  attitude-changes-first
  attitude-changes-second
  att-mean
  att-median
  att-skew
  positive-memos
  neutral-memos
  negative-memos
  plus-to-minus
  minus-to-plus
]

breed [ citizens citizen ]
breed [ channels channel ]
breed [ sources source ]
breed [ places place ]
breed [ dummies dummy ]

turtles-own
[
  num               ;; this is a counter for all turtles such that for each breed counting begins with 0
  attitude-position ;; this is the same as position on the issue as in AbeBer63, used in A3 and A10 for citizens
]                   ;; and attitude position unsed in A18 for sources

citizens-own
[
  memos-in-memory  ;; used in A1 ans A2
  memos-posted     ;; used in B24
  predisposition        ;; used in A9
  interest-in-the-issue ;; used in A2 and A12 and A13 and A14
  probability-to-vote
  vote                  ;; used in C1 and C2
  consistence-multiplier     ;; two auxiliar global variables used in
  predisposition-multiplier  ;; two procedures to avoid double calculation
  attitude-position-at-start
  colour-at-start
  attitude-change-first
  attitude-change-second
  changed-attitude?
  witness?                   ;; for test purposes
]

channels-own
[
  current-memos
  bias
]

places-own
[
  posted-memos
]

sources-own
[
  memos-given
]


directed-link-breed [ go-to goes-to ]                           ;; see p.113

directed-link-breed [ are-exposed-to is-exposed-to ]       ;; used in A1
are-exposed-to-own                ;; from citizen to source
[
  memo-match      ;; used in A5
  attitude        ;; used in A6
  receptivity     ;; used in A6
  satisfaction    ;; used in A15-A17 and A21 and A22
]

directed-link-breed [ are-used-by is-used-by ]             ;; used in A2
are-used-by-own                  ;; from channel to source
[
  attention-value
]

directed-link-breed [ are-attracted-by is-attracted-by ] ;; used in A1
are-attracted-by-own            ;; from citizen to channel
[
  attraction
]

directed-link-breed [ talk-to talks-to ]                         ;; used in B1
talk-to-own                      ;; from talking citizen to listening citizen
[
  memo-match ;; used in B5
  attitude
  receptivity ;; used in B2
]

to setup
  clear-all
  set attitude-changes-first []
  set attitude-changes-second []
  set att-mean 0.0
  set att-median 0.0
  set att-skew 0.0
  reset-timer
  reset-ticks
  set multi-run? behaviorspace-run-number > 0
;; setting names for output directories and output files
  set run-identifier format 4 behaviorspace-run-number
  let experiment-name ( word "ex-" n-sources "-" n-channels "-" n-places "-" n-citizens )
  if not pathdir:isDirectory? experiment-name [ pathdir:create experiment-name ]
  set directory-name ( word experiment-name "/run__" run-identifier )
;  set directory-name ( word experiment-name "/run_" rectify-date "_" run-identifier )
  if not pathdir:isDirectory? directory-name [ pathdir:create directory-name ]
  set logfile-name ( word directory-name "/log.txt" )
  set agentfile-name  ( word directory-name "/" run-identifier "-agents.csv")
  set sourcefile-name ( word directory-name "/" run-identifier "-sources.csv")
  if file-exists? logfile-name
  [
    file-close
    file-delete logfile-name
  ]
  file-open logfile-name
  file-print logfile-name
  file-print sourcefile-name
  file-close
  if file-exists? logfile-name [ print ( word " logfile is " logfile-name "." ) ]
  clear-output
;;  random-seed 47822               ;; only for test purposes
  output-print ( word run-identifier ":" n-channels "_" n-sources "_" timer ": setup starts " )
  set ERROR? false
  set test-level 20
  let ttl -1
  set max-memory-size-used 0
  set max-memory-size-current 0
  set min-memory-size-current 10000
  ask patches [ set pcolor white ]

  put-color-legend-to-view
  if n-sources < 1 or n-sources > 10 [ stop ]
  if n-channels < 1 or n-channels > 10 [ stop ]
  if n-citizens < 2 [ stop ]
  if n-places < 1 [ stop ]
  create-sources n-sources
  [
    set num who - count dummies
    set shape "coin heads"
    set color yellow
    setxy ( ifelse-value odd? n-sources [ ( 1 + num - floor ( n-sources / 2 ) ) * 3 ]
                                        [ ( 1 + num - floor ( n-sources / 2 ) ) * 3 - 1.5 ] ) ( max-pycor )
    set attitude-position random-normal 0.0 1.0
    set memos-given []
    convert-to-color attitude-position
  ]
  file-open sourcefile-name
  file-print "run. who, attitude"
  ifelse file-exists? sourcefile-name
  [ output-print ( word sourcefile-name " opened for source data." ) ]
  [
    output-print ( word sourcefile-name " not opened.")
    set ERROR? true
    stop
  ]
  ask sources [ file-print ( word n-channels " " n-sources " run_" run-identifier "," who " " num "," attitude-position  ) ]
  file-close
  create-channels n-channels
  [
    set num who - count dummies - n-sources
    set shape "computer workstation"
    set color blue
    setxy ( ifelse-value odd? n-channels [ ( 1 + num - floor ( n-channels / 2 ) ) * 3 ]
                                         [ ( 1 + num - floor ( n-channels / 2 ) ) * 3 - 1.5 ] ) ( max-pycor - 4 )
    set current-memos []
    create-are-used-by-to sources [ set attention-value ( 1 + random 2 ) / 2 ] ;; 0.5 for "low" 1 for "high" as in A4
  ]
  trace ttl ( word "channels to sources: " count are-used-by " links" )
  create-places n-places
  [
    set num who - count dummies - n-sources - n-channels
    set shape "flag"
    set color red
    move-to-position
    set posted-memos []
  ]
  create-citizens n-citizens
  [
    set num who - count dummies - n-sources - n-channels - n-places
    set shape "person"
    set color gray
    move-to-position
    create-go-to-to turtle-set min-n-of 2 places [ distance myself ]
    set interest-in-the-issue random-normal 0.0 1.0
    set attitude-position random-normal 0.0 1.0
    set attitude-position-at-start attitude-position
    convert-to-color attitude-position-at-start
    set colour-at-start color
    set attitude-change-first 0.0
    set attitude-change-second 0.0
    set changed-attitude? false
    set predisposition random-normal 0.0 1.0
    set memos-in-memory[]
    set memos-posted []
    set witness? false
    set vote vote-or-abstain
  ]
  file-open agentfile-name
  file-print "run_and_week, who, weight, predisposition, interest, attitude, att_change_1, att_change_2, prob_vote, vote"
  ifelse file-exists? agentfile-name
  [ output-print ( word run-identifier "_" n-channels "_" n-sources "_" timer " " agentfile-name " opened for agent data." ) ]
  [ output-print ( word "ERROR: " agentfile-name " could not be created " ) stop ]
  ask citizens
  [
    file-print ( word n-channels " " n-sources " run_" run-identifier "_week_--0," who ",1.0," predisposition
                      "," interest-in-the-issue "," attitude-position ","
      attitude-change-first "," attitude-change-second "," probability-to-vote "," vote )
  ]
  file-close
  trace ttl ( word "citzens to places: " count go-to " links" )
  ask go-to [ hide-link ]
  ask citizens
  [
    let potential-communicators turtle-set nobody
    ask out-goes-to-neighbors
    [
      let my-customers in-goes-to-neighbors
      set potential-communicators ( turtle-set potential-communicators my-customers )
    ]
    let n-of-potential-communicators count other potential-communicators
    if n-of-potential-communicators > 6 [ set n-of-potential-communicators 6 ]
    let communication-partners n-of n-of-potential-communicators other potential-communicators
    create-talk-to-to communication-partners
    create-are-attracted-by-to channels
    ask are-attracted-by
    [
      hide-link
      set attraction random-float 1.0
      convert-to-color attraction
    ]
  ]
  trace ttl ( word "citizens to citizens: " count talk-to " links" )
  trace ttl ( word "citizens to channels: " count are-attracted-by " links" )
  if test-level > 17 [ ask one-of citizens [ set witness? true ] ]
  set ERROR? false
  set att-mean mean [ attitude-position ] of citizens
  set att-median median [ attitude-position ] of citizens
  set att-skew skewness [ attitude-position ] of citizens
  setup-temporary-plot-pens
  update-plots
  output-print ( word run-identifier "_" n-channels "_" n-sources "_" timer ": setup complete " )
  reset-ticks
end

to go
  set attitude-changes-first []
  set positive-memos 0
  set neutral-memos 0
  set negative-memos 0
  if ERROR? [ show "ERROR" stop ]
  let ttl -1
  output-print ( word run-identifier "_" n-channels "_" n-sources "_" timer ": week " ticks " starts" )
  first-half-week
  if attitude-changes-first = 0 [ set attitude-changes-first [] ]
  second-half-week
  output-print ( word run-identifier "_" n-channels "_" n-sources "_" timer ": poll for week " ticks " done")
  set att-mean mean [ attitude-position ] of citizens
  set att-median median [ attitude-position ] of citizens
  set att-skew skewness [ attitude-position ] of citizens
  do-plot
  if not multi-run? [ export-interface ( word directory-name "/AbeBer63-" format 3 ticks ".png" ) ]
  tick
  if ticks > max-ticks [ stop ]
end ;; go

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; first half-week: citizens are exposed to sources via channels
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to first-half-week
  let ttl -1
  ask sources
  [
    let new-memo create-memo
    ask in-is-used-by-neighbors [ set current-memos accept-distributed-memos current-memos new-memo ]
  ]
  output-print ( word run-identifier "_" n-channels "_" n-sources "_" timer ": memos distributed" )
  ask channels [ trace-other ttl ( word "current memos " )  current-memos true ]
;; now channels are full of memos; citizens will now pick them up
  file-open logfile-name
  ask citizens
  [
    trace-citizen ttl ( word "memos-in-memory at begin of week: " length memos-in-memory " items in "  )
                         memos-in-memory true
    let this-citizen who
    let memos-exposed-to []
    let new-memos []
    let my-channel-links my-out-are-attracted-by
    trace-citizen ttl my-channel-links [] false
    let my-channels turtle-set remove-duplicates [ end2 ] of my-channel-links
    trace-citizen ttl my-channels [] false
    set memos-exposed-to get-memos-exposed-to my-channels new-memos
    trace-citizen ttl ( word "memos-exposed-to: " length memos-exposed-to " items in "  ) memos-exposed-to true
;; now citizens have listened to all memos from all channels;
;; they will now sort them out
    let memos-to-memory remove-memos-via-scorned-channels memos-exposed-to
    trace-citizen ttl ( word "memos-to-memory scorned removed, kept " length memos-to-memory " items in "  ) memos-to-memory true
;; rule A1 executed
    set memos-to-memory remove-memos-for-lack-of-interest memos-to-memory
    trace-citizen ttl ( word "memos-to-memory interesting kept: " length memos-to-memory " items in "  ) memos-to-memory true
;; rule A2 executed, broadcast-memo now contains all memos to which citizen was exposed in this time step
    let my-old-memos memos-in-memory
    set memos-in-memory remove-duplicates sentence memos-in-memory memos-to-memory
    trace-citizen ttl ( word "memos-in-memory duplicates removed, kept " length memos-in-memory " items in "  ) memos-in-memory true
;; memos-in-memory now contains all memos to which citizen was exposed from the beginning (except those forgotten or neglected)
;; receptivity of citizen with respect to each source has to be calculated
    set memos-in-memory calc-receptivity-of-memos memos-in-memory
;; then the memo-match has to be calculated for each is-exposed-to
    ask my-out-are-exposed-to [ set memo-match calc-memo-match end1 end2 ]
;; after applying rule A6 and A7 accepted memos have state +1 (and this state will never be changed),
;; not accepted memo have state -1, and this, too will never be changed
    trace-citizen ttl ( word "memos in memory after rule A6: "  )  memos-in-memory true
;; rules A3 to A10 executed
;; apply rule A11
    set memos-in-memory upd-memos-from-disliked-sources memos-in-memory
;; apply rules A12-A14
    set interest-in-the-issue changed-interest-in-the-issue my-out-are-exposed-to
;; apply rules A15-A20 about citizen's attitude and interest changes
;;;;;;;    apply-attitude-and-interest-change
    ask my-out-are-exposed-to [ set attitude changed-attitude-toward-source ]
    ask my-out-are-attracted-by [ set attraction changed-attitude-toward-channel ]
;; for each source to which this citizen is connected!
    let former-attitude attitude-position
    set attitude-position changed-own-position my-out-are-exposed-to
    set attitude-change-first attitude-position - former-attitude
    trace ttl ( word former-attitude " -> " attitude-position " ## " attitude-change-first )
    if former-attitude != attitude-position [ set attitude-changes-first fput ( attitude-change-first ) attitude-changes-first ]
    set changed-attitude? changed-attitude? or ( former-attitude != attitude-position )
;;
    set memos-in-memory remove-older-duplicates memos-in-memory
;;    convert-to-color attitude-position
  ]
  output-print ( word run-identifier "_" n-channels "_" n-sources "_" timer ": first half of week " ticks " over" )
  trace ttl ( word "mid-week memory sizes " sort [ length memos-in-memory ] of citizens )
end ;; first-half-week

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; second half-week: citizen are exposed to each other
;; citizen create memos and put them on the blackboards of the meeting places
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to second-half-week
  let ttl -1
  set max-memos-received 0
  ask places [ set posted-memos [] ]
  ask citizens
  [
    let new-memo []
    set memos-posted remove-older-duplicates memos-posted
    let reusable-memos filter [ m -> state m = 1 ] memos-posted
    trace-citizen ttl ( word length reusable-memos " reusable-memos "  ) reusable-memos true
    let ran3 random-float 1.0

    ifelse empty? reusable-memos
    [ set new-memo create-memo ]          ;; entirely new memo
    [
      ifelse ran3 < 0.9
      [
        trace-citizen ttl ( word ran3 ": " length reusable-memos " reusable-memos "  ) reusable-memos true

        set new-memo create-memo ]        ;; entirely new memo with p=0.5
      [                                             ;; or
        trace-citizen ttl ( word ran3 ": " length reusable-memos " reusable-memos "  ) reusable-memos true

        set new-memo one-of reusable-memos
        set new-memo set-state new-memo 0 ;; reused memo, state reset to pending (0)
      ]
    ]
    trace-citizen ttl ( word "current reusable memos "  ", new or reused " new-memo ) reusable-memos true

    set memos-posted lput new-memo memos-posted
    ask out-goes-to-neighbors [ set posted-memos accept-posted-memos posted-memos new-memo ]
  ]
  output-print ( word run-identifier "_" n-channels "_" n-sources "_" timer ": citizens have posted" )
  ask citizens
  [
    let this-citizen who
    let memos-seen []
    let my-place-links my-out-go-to
    trace-citizen ttl ( word "current memory size (1) " length memos-in-memory ) [] false
    trace-citizen ttl my-place-links [] false
    let my-places turtle-set remove-duplicates [ end2 ] of my-place-links
    trace-citizen ttl my-places [] false
    ask my-places
    [
      set memos-seen sentence posted-memos memos-seen
      trace ttl ( word "read " length posted-memos ", total " length memos-seen )
    ]
    trace-citizen ttl ( word "memos seen: " length memos-seen )  memos-seen true
    set memos-in-memory sentence remove-duplicates memos-in-memory memos-seen
    trace-citizen ttl ( word "current memory size (2) " length memos-in-memory ) [] false
    set memos-in-memory remove-memos-for-lack-of-interest memos-in-memory
    trace-citizen ttl ( word "current memory size (3) " length memos-in-memory ) [] false
    set memos-in-memory remove-older-duplicates memos-in-memory
    trace-citizen ttl ( word "current memory size (4) " length memos-in-memory ) [] false
;; rule B1 executed, memos-received now contains all posted memos which citizen has read
    ask my-out-talk-to [ set memo-match calc-memo-match end1 end2 ] ;; B5
    set memos-in-memory calc-receptivity-of-memos memos-in-memory
    trace-citizen ttl ( word "memos in memory after rule 6: " length memos-in-memory ) [] false
;; rules A3 to A10 executed
;; apply rule 11
    set memos-in-memory upd-memos-from-disliked-sources memos-in-memory
    trace-citizen ttl ( word "memos in memory after rule 11: " length memos-in-memory ) [] false
;; apply rules A12-A14
    set interest-in-the-issue changed-interest-in-the-issue my-out-talk-to

    trace-citizen ttl ( word "memos in memory after rule 14: " length memos-in-memory ) [] false
  ]
  trace ttl ( word "memory after exposure to peers " sort [ length memos-in-memory ] of citizens )
  output-print ( word run-identifier "_" n-channels "_" n-sources "_" timer ": citizens read postings" )
  ask citizens
  [
    let this-citizen self
    let former-attitude attitude-position
;; apply rules A15-A20 about citizen's attitude and interest changes
    ask my-out-are-exposed-to
    [
      set attitude changed-attitude-toward-source
      set satisfaction changed-satisfaction-with-source
    ]
    convert-to-color attitude-position
;; rules B6-B21 applied
;; B22-B24 about forgetting memos to follow
    if not empty? memos-in-memory
    [
      trace-citizen ttl ( word length memos-in-memory " before removing duplicates" ) [] false
      set memos-in-memory remove-older-duplicates memos-in-memory
      trace-citizen ttl ( word length memos-in-memory " after removing duplicates, before forgetting" ) [] false
      set memos-in-memory forget-memos memos-in-memory
      trace-citizen ttl ( word length memos-in-memory " after forgetting" ) [] false
    ]
  ]
  output-print ( word run-identifier "_" n-channels "_" n-sources "_" timer ": citizens updated attitudes towards sources and memories" )
  trace ttl ( word "memory sizes after forgetting " sort [ length memos-in-memory ] of citizens )
;; B25-27 about position and interest changes
  ask citizens
  [
    let lar length memos-in-memory
    if lar > max-memos-received [ set max-memos-received lar ]
    trace-citizen ttl ( word "memos-received " lar ", max " max-memos-received ) [] false
  ]
  ask channels
  [
    set bias ( mean map [ m -> opinion m ] current-memos + 1 ) / 2
  ]
  set attitude-changes-second []
  let max-memos-in-memory max [ length memos-in-memory ] of citizens
  ask citizens with [ not empty? memos-in-memory ]
  [
;; B25
    let rank ( length memos-in-memory ) / max-memos-received
    let long-opinion-list map [ m -> opinion m ] memos-in-memory
    let short-opinion-list filter [ m -> state m = 1 ] memos-in-memory
    let opinion-list map [ m -> opinion m ] filter [ m -> state m = 1 and opinion m != 0 ] memos-in-memory
    let old-attitude attitude-position
    if not empty? opinion-list
    [
      let majority-memo-positions sum opinion-list
      set majority-memo-positions majority-memo-positions / max-memos-in-memory
;;      if majority-memo-positions > 0 [ set majority-memo-positions 1 / length opinion-list ]
;;      if majority-memo-positions = 0 [ set majority-memo-positions 0 ]
      set attitude-position upd-attitude attitude-position majority-memo-positions rank
      set attitude-change-second attitude-position - old-attitude
      if attitude-change-second != 0 [ set attitude-changes-second fput ( attitude-change-second ) attitude-changes-second ]
      if abs attitude-change-second > 0.1 [ trace ttl ( word "abs " attitude-change-second ">0.1!!!" ) ]
      set changed-attitude? changed-attitude? or ( old-attitude != attitude-position )
      trace-citizen ttl ( word " B25: " long-opinion-list opinion-list " " sum opinion-list " mmim " max-memos-in-memory ) short-opinion-list true
      trace-citizen ttl ( word " B25: attitude changed from " old-attitude " towards " majority-memo-positions
           " by " attitude-change-second " to " attitude-position ", mmp " majority-memo-positions ", rank " rank ) [] false
      trace ttl "==========================="
    ]
;; B26
    set interest-in-the-issue upd-interest interest-in-the-issue rank
;; B27
    ask my-out-are-attracted-by [ set attraction upd-attraction attraction [ attitude-change-second ] of myself ]
  ]
  output-print ( word run-identifier "_" n-channels "_" n-sources "_" timer ": citizens updated positions" )
;; C1-C2 to follow
  file-open agentfile-name
  output-print ( word run-identifier " " n-channels " " n-sources " " agentfile-name " opened for agent data." )
  ask citizens
  [
    set vote vote-or-abstain
    file-print ( word n-channels ":" n-sources "_run_" run-identifier "_week_" format 3 ticks "," who ",1.0," predisposition ","
      interest-in-the-issue "," attitude-position "," attitude-change-first "," attitude-change-second "," probability-to-vote "," vote )
  ]
  file-flush
  file-close
end ;; second-half-week

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; function CuS -> A creating a new memo for either a citizen
;; or a source (the caller of this function)
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to-report create-memo
  let ttl -1
  let this-opinion 0
  let this-aspect item ( random 5 ) "abcde"
  ( ifelse
    attitude-position < -0.430723756906 [ set this-opinion -1 set negative-memos negative-memos + 1 ]
    attitude-position < 0.430723756906  [ set this-opinion 0 set neutral-memos neutral-memos + 1 ]
    [ set this-opinion 1 set positive-memos positive-memos + 1 ]
  )
  trace ttl ( word "attitude-position: " attitude-position " opinion " this-opinion )

  let an-memo ( list "from " who               ;; the memo is a sextuple consisting of source-id,
                     " via " 0            ;; the channel-id through which it came
                     " opi " this-opinion ;; the opinion being pro (+1), con (-1) or neutral (0)
                     " asp " this-aspect  ;; the aspect of the fluoridation issue currently mentioned, one of the first five letters of the alphabet
                     " at " ticks         ;; the time when it was sent
                     " sta " 0            ;; the state (0: pending, 1 accepted by the citizenwhich has this memo in its memory, -1 rejected
                     " fgt " 0 )          ;; used only in B22
  trace ttl ( word "my channels: " in-is-used-by-neighbors )
  report an-memo
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; function A^n x A -> A^(n+1)
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to-report accept-distributed-memos [ old-memos new-memo ]
  let ttl -1
  let old-memo new-memo
  set new-memo set-broadcaster new-memo who
  trace ttl ( word "channel " who " changed memo from " old-memo " to " new-memo )
  report lput new-memo old-memos
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; function A^n x A -> A^(n+1)
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to-report accept-posted-memos [ old-memos new-memo ]
  let ttl -1
  trace ttl "entering post-memos"
  trace ttl ( word "my places are " [ who ] of out-goes-to-neighbors )
  let old-memo new-memo
  set new-memo set-broadcaster new-memo -1
  trace ttl ( word "place " who " changed memo from " old-memo " to " new-memo )
  report remove-duplicates lput new-memo old-memos
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; function ø->A^n getting memos from all channels
;; mentioned in A1 and A2
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to-report get-memos-exposed-to [ these-channels new-memos ]
  ask these-channels
  [
    set new-memos remove-duplicates sentence current-memos new-memos
  ]
  report new-memos
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; function A^n->A^m (m <= n) removing older similar memos from the same source
;; TODO
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to-report remove-older-duplicates [ memos ]
  if empty? memos [ report memos ]
  let ttl -1
  trace-citizen ttl ( word "looking for older duplicates in "  " (" length memos ")" ) memos true
  set memos sort-by
  [ [ element next-element ] -> reduce [ [ a b ] -> ( word a b ) ]  element < reduce [ [ c d ] -> ( word c d ) ] next-element
  ] memos
  trace-citizen ttl ( word "sorted memo list "  " (" length memos ")" ) memos true
;; now, the older of two adjacent similar memos needs to be removed
  let n-memos length memos
  let positions n-values ( n-memos - 1 ) [ i -> i ]
;;  file-open logfile-name
;;  file-show positions
  foreach positions
  [ i ->  ifelse is-similar? item i memos item ( i + 1 ) memos
    [
;;      file-show ( word i ": " item i memos ", " ( i + 1 ) ": " item ( i + 1 ) memos ": similar" )
      set memos replace-item i memos item ( i + 1 ) memos
    ]
    [
;;      file-show ( word i ": " item i memos ", " ( i + 1 ) ": " item ( i + 1 ) memos ": different" )
    ]
  ]
  trace-citizen ttl ( word "memo list replaced duplicates"  " (" length memos ")" ) memos true
  set memos remove-duplicates memos
  trace-citizen ttl ( word "memo list without older duplicates"  " (" length memos ")" ) memos true
  report memos
end

to-report is-similar? [ test-memo tested-memo ]
  let result false
  set result sender test-memo = sender tested-memo
    and broadcaster test-memo = broadcaster tested-memo
    and opinion test-memo = opinion tested-memo
    and aspect test-memo = aspect tested-memo
    and state test-memo = state tested-memo
    and week test-memo < week tested-memo
    and forgettability tested-memo >= 0
    and forgettability test-memo >= 0
  report result
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; function A^n->A^m (m <= n) removing memos received via unattractive channels
;; this function realises rule A1
;; TODO
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to-report remove-memos-via-scorned-channels [ memos ]
  let ttl -1
  let memos-from-attractive-channels []
  foreach memos [ this-memo -> set memos-from-attractive-channels remove-scorned this-memo memos-from-attractive-channels ]
  trace-citizen ttl ( word length memos-from-attractive-channels " out of " length memos " received from attractive channels " )  [] false
  report memos-from-attractive-channels
end

to-report remove-scorned [ this-memo mfac ]
  let ttl -1
  let this-channel broadcaster this-memo
  let this-channel-link is-attracted-by who this-channel
  ifelse this-channel-link != nobody
  [
    let this-channel-attractivity [ attraction ] of this-channel-link
    let ran random-float 1.0
    if this-channel-attractivity < ran
    [
      trace-citizen ttl ( word "keeping " this-memo " for sufficient attractivity of channel: " this-channel-attractivity " < " ran )  [] false
      set mfac fput this-memo mfac
    ]
  ]
  [
    trace-citizen ttl ( word "removing " this-memo " for missing link to channel" )  [] false
  ]
  report mfac
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; function A^n->A^m (m <= n) removing memos received for lack of interest
;; this function realises rule A2
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to-report remove-memos-for-lack-of-interest [ memos ]
  let interesting-memos []
  let ttl -1
  let lba length memos
  if lba = 0 [ report memos ]
  foreach memos [ this-memo -> set interesting-memos find-interesting-memos interesting-memos this-memo ]
  let lba2 length interesting-memos
  trace-citizen ttl ( word "received " lba " memos at " ticks " with interest " interest-in-the-issue ", kept " lba2 " = " ( lba2 / lba )  )  [] false
  report interesting-memos
end

to-report find-interesting-memos [ interesting-memos this-memo ]
  let ttl -1
  let ran random-normal 0.0 1.0
  ifelse ran > interest-in-the-issue ;; i.e. with probability interest-in-the-issue this memo is moved to interesting-memos
  [
    trace-citizen ttl ( word "memo " this-memo " is sufficiently interesting: " ran ">" interest-in-the-issue )  [] false
    set interesting-memos fput this-memo interesting-memos
  ]
  [
    trace-citizen ttl (word "memo " this-memo " is not sufficiently interesting: " ran "<=" interest-in-the-issue ) [] false
  ]
  let this-citizen who
  let this-sender sender this-memo
  ifelse broadcaster this-memo >= 0 ;; i.e. sender is a source, not a citizen
  [ create-is-exposed-to-to source ( sender this-memo ) [ set attitude random-normal 0.0 1.0 ] ]
  ;; else: sender is another citizen, but not this one
  [
    if this-citizen != this-sender ;; i.e. this-sender is another citizen
    [
      let this-link talks-to this-sender this-citizen
      trace-citizen ttl ( word "talks-to from " this-sender " to " this-citizen " becomes active, " this-link )  [] false
      if this-link != nobody [ ask this-link [ if color = gray [ set color white ] ] ]
    ]
  ]
  report interesting-memos
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; function A^n -> A^n calculates the receptivity of memos and resets their state element accordingly
;; realises rules A3, A6, A7, A8, A9, A10
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to-report calc-receptivity-of-memos [ memos ]
  let ttl -1
  file-open logfile-name
  if witness? [ document-memos memos ]
  foreach memos [ this-memo -> set memos calc-rec-memos this-memo memos ]
  if witness? [ document-memos memos ]
  report memos
end

to-report calc-rec-memos [ this-memo memos ]
  let ttl -1
  let this-link nobody
  set memos remove this-memo memos
  let this-source sender this-memo
  let attitude-towards-this-source 0.0
  ifelse broadcaster this-memo < 0
  [ set this-link talks-to this-source who ]
  [ set this-link is-exposed-to who this-source ]
  if this-link != nobody [ ask this-link [ set attitude-towards-this-source attitude ] ]
  let this-channel broadcaster this-memo
  let this-opinion opinion this-memo
  let this-aspect aspect this-memo
  let this-state state this-memo
  let this-source-position 0
  let this-source-attention-value 1
  ifelse this-channel >= 0
  [
    set this-source-position [ attitude-position ] of source this-source
    set this-source-attention-value [ attention-value ] of is-used-by this-channel this-source ;; A4
  ]
  [ set this-source-position [ attitude-position ] of citizen this-source ];; finding the attention-value relevant for the receptivity of this-memo
;; calculating the extremity of citizen's attitude position
  let extremity abs attitude-position / 4   ;; A3
;; calculating the receptivity of this memo from extremity, attention-value and memo-match
  let this-receptivity 0
  if this-link != nobody
  [
    set this-receptivity calc-receptivity extremity this-source-attention-value [ memo-match ] of this-link
    ask this-link [ set receptivity this-receptivity ]
  ]
;; memos are now analysed for being accepted or not, except their state was already +1 according to A7
;; the probability of being accepted is calculated from the citizen's attitude toward the source and its receptivity to the source (A6)
  let probability-of-being-accepted this-receptivity
  if state this-memo != 1
  [
;; apply rule 8
    let similar-memos find-similar-memos ( opinion this-memo ) ( aspect this-memo ) memos ;; rule A8
    let ran2 random-float 1.0
    if similar-memos = nobody ;; rule A8: this-memo is entirely new wrt aspect and opinion -> receptivity higher
    [
      set probability-of-being-accepted increase probability-of-being-accepted 0.2
    ]
    if similar-memos != nobody ;; rule A8: there were earlier negative memos wrt this aspect -> receptivity lower
    [
      set probability-of-being-accepted decrease probability-of-being-accepted 0.2
    ]
;; apply rule 9
    if abs ( opinion this-memo - predisposition ) > 0.5 [ set probability-of-being-accepted decrease probability-of-being-accepted 0.2 ]
;; apply rule 10
    trace ttl ( word "A10: state of this memo is " state this-memo ", " this-memo )
    ifelse state this-memo != 0
    [
      trace-citizen ttl ( word "A10: no state change as state was " state this-memo ) [] false
    ]
    [
      if abs ( opinion this-memo - attitude-position ) > 0.5 [ set probability-of-being-accepted decrease probability-of-being-accepted 0.2 ]
      ifelse probability-of-being-accepted > ran2  ;; rule A6: for, e.g., ran2=0.8 the probability of setting state=-1 is 0.2
      [
        set this-memo set-state this-memo -1
        trace-citizen ttl ( word  "A10: " this-memo " changed from " this-state " to -1 because of " probability-of-being-accepted " > " ran2 )  [] false
      ]
      [
        set this-memo set-state this-memo 1
        trace-citizen ttl ( word  "A10: " this-memo " changed from " this-state " to +1 because of " probability-of-being-accepted " <= " ran2 )  [] false
      ]
      trace-citizen ttl ( word "A10: overall change of probability-of-being-accepted from " this-receptivity " to " probability-of-being-accepted )  [] false
    ]
  ]
  set memos fput this-memo memos
  report memos
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; function A^n -> A^n changing the state of memos which came from disliked sources
;; implements A11
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to-report upd-memos-from-disliked-sources [ memos ]
  let ttl -1
  let disliked-sources-memos find-disliked-sources-memos memos ;; rule 11
  let changed-memos []
  trace-citizen ttl ( word "memos from disliked sources are " length disliked-sources-memos )  [] false
  if disliked-sources-memos != nobody
  [
    foreach memos [ test-memo -> set changed-memos upd-mem-disliked test-memo disliked-sources-memos changed-memos ]
  ]
  report changed-memos
end

to-report upd-mem-disliked [ test-memo disliked-sources-memos changed-memos ]
  foreach disliked-sources-memos [ tested-memo -> set changed-memos upd-mem-disliked-tested tested-memo test-memo changed-memos ]
  report changed-memos
end

to-report upd-mem-disliked-tested [ tested-memo test-memo changed-memos ]
  let ttl -1
;; tested-memo ?       ;; came from disliked source
  trace-citizen ttl ( word "memos to compare are " test-memo " and " tested-memo )  [] false
  if ( aspect test-memo ) = ( aspect tested-memo ) ;; same aspect
  and ( opinion test-memo ) = (opinion tested-memo ) ;; same opinion
  and ( sender test-memo ) != ( sender tested-memo ) ;; different source
  [
    set tested-memo set-state tested-memo -1 ;; the disliked-source has the same opinion about this aspect
                                             ;; as the source of this tested-memo, hence this citizen
                                             ;; rejects this memo
    trace-citizen ttl ( word "set state to -1 for " test-memo )  [] false
  ]
  set changed-memos remove-duplicates fput test-memo changed-memos
  report changed-memos
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; function C x S -> R changes the interest in the issue because of memo match
;; implements A12--A14
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to-report changed-interest-in-the-issue [ these-links ]
  let ttl -1
  let n-of-received-memos length memos-in-memory
  let new-interest interest-in-the-issue
;; A12-A22 to follow
  ask these-links
  [
    let this-source [ who ] of end2
;; rule A12
    if memo-match > 0.5
    [
      let this-memo-match memo-match
      ask end1
      [
        let old-interest interest-in-the-issue
        set interest-in-the-issue increase interest-in-the-issue ( this-memo-match - 0.5 )
        ifelse old-interest != interest-in-the-issue
        [ trace-citizen ttl ( word "A12 changed interest because of source " this-source " from " old-interest " to " interest-in-the-issue )  [] false ]
        [ trace-citizen ttl ( word "A12: no change of interest because of source " this-source " as memo match is " this-memo-match )  [] false ]
      ]
    ]
    if n-of-received-memos > 0
;; rule A13
    [
      ask end1
      [
        let n-of-relevant-memos length filter [ m -> sender m = this-source and state m = 1 ] memos-in-memory
        set predisposition-multiplier  n-of-relevant-memos / n-of-received-memos
        let old-interest interest-in-the-issue
        set new-interest interest-in-the-issue
        set new-interest increase new-interest predisposition-multiplier
        ifelse old-interest != new-interest
        [ trace-citizen ttl ( word "A13 changed interest because of source " this-source " from " old-interest " to " new-interest )  [] false ]
        [ trace-citizen ttl ( word "A13: no change of interest because of source " this-source " as predisposition m. is " predisposition-multiplier )  [] false ]
;; rule A14
        let n-of-inconsistent-memos count-inconsistent-memos memos-in-memory
        set consistence-multiplier n-of-inconsistent-memos / n-of-received-memos / 10
        set old-interest interest-in-the-issue
        set new-interest decrease new-interest consistence-multiplier
        ifelse old-interest != new-interest
        [ trace-citizen ttl ( word "A14 changed interest because of source " this-source " from " old-interest " to " new-interest ) [] false ]
        [ trace-citizen ttl ( word "A14: no change of interest because of source " this-source " as consistence m. is " consistence-multiplier )  [] false ]
      ]
    ]
  ]
  report new-interest
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; function C x S -> R changes the attitude of citizen towards source
;; implements A21
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to-report changed-attitude-toward-source
  let ttl -1
  trace ttl "changed-attitude-toward-source starts"
  let old-attitude attitude
;; A21 about change of attitude towards source
  ifelse satisfaction > 0.5
  [ set attitude increase attitude ( satisfaction * 0.5 ) ]
  [ set attitude decrease attitude ( satisfaction * 0.5 ) ]
  trace ttl ( word end1 "'s attitude towards " end2 " was " old-attitude ", is now " attitude )
  report attitude
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; function C x S -> R changes a citizen's own position on the issue
;; implements A18--A20
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to-report changed-own-position [ these-links ]
  let ttl -1
;; A18: citizen's attitude changes only if source's attitude is more extreme than the citizen's
  let attitude-changes-possible? false
  let new-position attitude-position
;; rule A18
  ask these-links
  [
    let my-attitude [ attitude-position ] of end1
    let its-attitude [ attitude-position ] of end2
    let this-consistence-multiplier [ consistence-multiplier ] of end1
    let this-predisposition-multiplier [ predisposition-multiplier ] of end1
    let my-memos-received [ memos-in-memory ] of end1
    set attitude-changes-possible? more-extreme? my-attitude its-attitude
    trace ttl ( word "A18-A20 changed-own-position, my-attitude " my-attitude ", its attitude " its-attitude ", " attitude-changes-possible? )
    let memos-from-this-source filter [ m -> sender m = [ who ] of end2 ] my-memos-received
;; A18-A20 about changes of attitude-position
    if attitude-changes-possible?
    [
      let attitude-position-change satisfaction ;; A19 postulates a "direct function"
      set new-position 0
      let value ( abs ( my-attitude - its-attitude ) * attitude-position-change )
      ifelse attitude < 0.0
      [ set new-position increase my-attitude value ]
      [ set new-position decrease my-attitude value ]
      trace ttl ( word "A18-A20 attitude toward source " attitude ", old position " my-attitude ". source's position " its-attitude
         ", change factor (satisfaction) " attitude-position-change ", new position " new-position )
    ]
    let attitude-change new-position - my-attitude
    trace ttl ( word "A18-A20: " [ who ] of end1 "'s attitude was " my-attitude ", is now " new-position ", its-attitude: " its-attitude
      ifelse-value attitude-change > 0 [ ( word ", change was: " attitude-change ) ] [ " no change " ] )
  ]
  report new-position
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; function C x S -> (0,1) changes the satisfaction of citizen with respect to a source
;; implements A15--A17
;; TODO, needs check
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to-report changed-satisfaction-with-source
  let this-consistence-multiplier [ consistence-multiplier ] of end1
  let this-predisposition-multiplier [ predisposition-multiplier ] of end1
;; A15-A17 about satisfaction with exposure
  let new-satisfaction memo-match ;; A15 postulates a "direct function"
  let my-attitude [ attitude-position ] of end1
  let memos-from-this-source filter [ m -> sender m = [ who ] of end2 ] [ memos-in-memory ] of end1
;; the next block needs test
;;  show ( word "nsf before: " new-satisfaction " tcm: " this-consistence-multiplier " tpm: " this-predisposition-multiplier )
  foreach memos-from-this-source
  [ ? -> set new-satisfaction change-satisfaction ? this-consistence-multiplier this-predisposition-multiplier my-attitude new-satisfaction ]
;;  show ( word "nsf after : " new-satisfaction )
  report new-satisfaction
end

to-report change-satisfaction [ m tcm tpm my-attitude new-satisfaction ]
  let consistence abs ( my-attitude - opinion m ) ;; A17
  ifelse consistence < 0.5
  [ set new-satisfaction increase new-satisfaction tcm ]
  [ set new-satisfaction decrease new-satisfaction tcm]
  set new-satisfaction decrease new-satisfaction tpm ;; A16
;;  show ( word "c: " consistence " nsf: " new-satisfaction )
  report new-satisfaction
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; function C x X -> R
;; A22 about change of attitude towards channel
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to-report changed-attitude-toward-channel
  let this-citizen [ who ] of end1
  let my-memos-received [ memos-in-memory ] of end1
  let this-channel [ who ] of end2
  let memos-through-this-channel filter [ m -> broadcaster m = this-channel ] my-memos-received
  let channel-satisfaction 0
  foreach memos-through-this-channel [ m ->  set channel-satisfaction changed-att-tow-chan m this-citizen channel-satisfaction ]
  report channel-satisfaction
end

to-report changed-att-tow-chan [ m this-citizen channel-satisfaction ]
  let this-satisfaction 0
  ask is-exposed-to this-citizen ( sender m ) [ set this-satisfaction satisfaction ]
  if this-satisfaction > channel-satisfaction [ set channel-satisfaction  this-satisfaction ]
  convert-to-color channel-satisfaction
  report channel-satisfaction
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; function A^n->A^m (m <= n) removing memos which can be forgotten for various reasons
;; implements B22
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to-report forget-memos [ memos ]
  let ttl -1
  if length memos < 5 [ report memos ]
  let this-citizen who
  let memos-to-be-kept []
  let memos-to-be-finally-kept []
  foreach memos [ tested-memo -> set memos-to-be-kept keep-memos tested-memo memos-to-be-kept memos ]

  set memos-to-be-kept remove-duplicates memos-to-be-kept
  let forgettable-memos filter [ m -> forgettability m > 0 ] memos-to-be-kept
  let n-forgettable-memos length forgettable-memos
  trace-citizen ttl ( word "end: " n-forgettable-memos " memos out of " length memos " / " length memos-to-be-kept " forgettable " ) [] false

  foreach memos-to-be-kept [ this-memo -> set memos-to-be-finally-kept keep-memos-finally this-memo memos-to-be-finally-kept memos-to-be-kept ]
  set memos-to-be-finally-kept remove-duplicates memos-to-be-finally-kept
  trace-citizen ttl ( word length memos-to-be-finally-kept " after forgetting (raw)" )  [] false
  set memos-to-be-finally-kept filter [ m -> forgettability m > 0 ] memos-to-be-finally-kept
  trace-citizen ttl ( word length memos-to-be-finally-kept " after forgetting (final)" )   [] false
  let latbfk length memos-to-be-finally-kept
  if latbfk > max-memory-size-used [ set max-memory-size-used latbfk ]
  if latbfk > 20
  [
    trace-citizen ttl ( word "forgetting all but 20 out of " length memos-to-be-finally-kept " memos " )  [] false
    set memos-to-be-finally-kept n-of 20 memos-to-be-finally-kept
  ]
  set latbfk length memos-to-be-finally-kept
  if latbfk > max-memory-size-current [ set max-memory-size-current latbfk ]
  if latbfk < min-memory-size-current [ set min-memory-size-current latbfk ]
  report memos-to-be-finally-kept
end

to-report keep-memos-finally [ this-memo memos-to-be-finally-kept memos-to-be-kept ]
  let ttl -1
  let position-of-this-memo position this-memo memos-to-be-kept
  let this-forgettability forgettability this-memo
  let ran random-float 1.0
  trace-citizen ttl ( word this-memo " could be forgotten if " this-forgettability " > " ran )  [] false
  if ran < this-forgettability
  [
    set this-memo set-forgettability this-memo -1
    trace-citizen ttl ( word "     forgets " this-memo " at " position-of-this-memo ", forgettability " this-forgettability " > " ran )  [] false
  ]
  set memos-to-be-finally-kept fput this-memo memos-to-be-finally-kept
  report memos-to-be-finally-kept
end


to-report keep-memos [ tested-memo memos-to-be-kept memos ]
  let ttl -1
  let this-citizen who
  let original-memo tested-memo
  let position-of-tested-memo position tested-memo memos
  let this-forgettability forgettability tested-memo
  let new-forgettability 0
  let this-state state tested-memo
  ifelse this-forgettability > 0
  [ set new-forgettability increase this-forgettability 0.75 ]
  [ set new-forgettability  0.75 ]
  set tested-memo set-forgettability tested-memo new-forgettability
;; B22
  if state tested-memo = 1 and attitude-position < 0.5 [ set tested-memo set-forgettability tested-memo 0.95 ]
  if state tested-memo = -1 and attitude-position > 0.5 [ set tested-memo set-forgettability tested-memo 0.95 ] ;; B22
  trace-citizen ttl ( word "forgettability changed from " this-forgettability " to " forgettability tested-memo " for " this-state ", ap " attitude-position )  [] false
;; B23
  if state tested-memo = 1 and predisposition < 0.5 [ set tested-memo set-forgettability tested-memo 0.95 ]
  if state tested-memo = -1 and predisposition > 0.5 [ set tested-memo set-forgettability tested-memo 0.95 ]   ;; B23

  trace-citizen ttl ( word "forgettability changed from " this-forgettability " to " forgettability tested-memo " for " this-state ", pr " predisposition )  [] false
  if state tested-memo = 1 ;;; questionable !!
  [
    let this-source sender original-memo
    if is-citizen? this-source
    [
      ask citizen this-source
      [
        trace-citizen ttl ( word "memos before change: "  ) memos-in-memory true
        if member? original-memo memos-in-memory
        [
          let its-memo-position position original-memo memos-in-memory
          let updated-memo original-memo
          trace-citizen ttl ( word this-source " will not forget " original-memo " because " this-citizen " had accepted it" ) [] false
          set updated-memo set-forgettability updated-memo 0.25
          set memos-in-memory replace-item its-memo-position memos-in-memory updated-memo
          trace-citizen ttl ( word "memo " its-memo-position " in " length memos " unforgettable" )  [] false
        ]
        trace-citizen ttl ( word "memos after  change: " memos-in-memory ) memos-in-memory true
      ]
    ]
  ]
  set memos-to-be-kept fput tested-memo memos-to-be-kept
  report memos-to-be-kept
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; function C-> R changing interest-in-the-issue according to B26
;; implements B26
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to-report upd-interest [ interest rank ]
  let ttl -1
  let old-interest interest
  let new-interest interest
  ifelse rank < 0.5
  [ set new-interest decrease new-interest ( ( 1 - rank ) / 4 ) ]
  [ set new-interest increase new-interest rank ]
  trace-citizen ttl ( word "B26 changed interest from " old-interest " to " new-interest ", rank was " rank) [] false
  report new-interest
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; function C-> R changing interest-in-the-issue according to B27
;; implements B27
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to-report upd-attraction [ attr this-change ]
  let att attr
  ifelse this-change > 0
  [ set att increase att ( this-change + 1 ) / 2 ]
  [ set att decrease att ( this-change + 1 ) / 2 ]
  report att
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; function CxS-> (0,1) calculating the memo match between a citizen and a source
;; implements A5
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to-report calc-memo-match [ this-citizen this-source ]
  let ttl -1
  let memos-from-this-source filter [ m -> sender m = [ who ] of this-source ] [ memos-in-memory ] of this-citizen
  let n-of-memos-from-this-source length memos-from-this-source
  let agreements 0
  let disagreements 0
  trace-other ttl ( word "all memos " length [ memos-in-memory ] of this-citizen " " ) [ memos-in-memory ] of this-citizen true
  trace-other ttl ( word "memos from " this-source ": " n-of-memos-from-this-source " " ) memos-from-this-source  true
  set agreements    length filter [ m -> state m > 0 ] memos-from-this-source
  set disagreements length filter [ m -> state m < 0 ] memos-from-this-source
  let memo-match-raw agreements - disagreements
  trace ttl ( word "memo-match raw: " memo-match-raw " agreements " agreements ", disagreements " disagreements )
  let result 0.0
  if n-of-memos-from-this-source  > 0
  [ set result memo-match-raw / n-of-memos-from-this-source ]
  trace ttl( word "memo match " result " agreements: " agreements ", disagreements: " disagreements ", memos: " n-of-memos-from-this-source )
  convert-to-color result
  report result
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; function calculating the receptivity, yields a value between 0 and 1
;; and is the higher the less extreme cpos is, the higher attval is
;; and the higher assmatch is
;; implements part of A3--A10, is called by calc-receptivity-of-memos
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to-report calc-receptivity [ extr attval assmatch ] ;; yields a value between 0 and 1 and is the higher
                                                     ;; the less extreme cpos is, the higher attval is
                                                     ;; and the higher assmatch is
  let ttl -1
  let return-value ( 1 - extr ) * attval * assmatch
  trace ttl ( word "extr " extr " attval " attval " assmatch " assmatch " receptivity " return-value )
  report return-value
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; function returning a sublist of memos with the specified opinion and aspect values
;; implements part of A8, is called by calc-receptivity-of-memos
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to-report find-similar-memos [ an-opinion an-aspect memos ]
  let similar-memos filter [ m -> opinion m = an-opinion and aspect m = an-aspect and state m = -1 ] memos
  report similar-memos
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; function A^n -> A^m, m <= n, returns memos from disliked sources
;; implements part of A11
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to-report find-disliked-sources-memos [ memos ]
  let ttl -1
  let disliked-memos []
  foreach memos [ this-memo -> set disliked-memos find-disl-sou-memos this-memo memos disliked-memos ]
  trace-citizen ttl ( word length disliked-memos " memos marked as disliked from "  length memos )  [] false
  report disliked-memos
end

to-report find-disl-sou-memos [ this-memo memos disliked-memos ]
  let ttl -1
  let this-source sender this-memo
  let this-channel broadcaster this-memo
  let this-link nobody
  ifelse this-channel < 0
  [ set this-link talks-to who this-source ]
  [ set this-link is-exposed-to who this-source ]
  if this-link != nobody
  [
    trace-citizen ttl ( word this-link " attitude " [ attitude ] of this-link ", receptivity " [ receptivity ] of this-link )  [] false
    ifelse [ attitude ] of this-link  < 0.5 and [ receptivity ] of this-link < 0.5
    [
      trace-citizen ttl ( word "memo " this-memo " marked as from disliked source or citizen " this-source )  [] false
      set disliked-memos fput this-memo disliked-memos
    ]
    [
      trace-citizen ttl ( word "memo " this-memo " marked as liked source or citizen " this-source )  [] false
    ]
  ]
  report  disliked-memos
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; function CxA->N counting memos in a citizen's memory which are inconsistent with its opinion
;; implements A11
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to-report count-inconsistent-memos [ memos ]
  report length filter [ m -> abs ( attitude-position - opinion m ) > 0.5 ] memos
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; function RxR-> boolean calculating whether the second argument is more extrem
;; than the first argument,reports an error if the first argument is out of bounds
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to-report more-extreme? [ inner outer ] ;;??? to be checked for sign equality in A18-A20
  let return-value? false
  if abs ( inner ) < abs ( outer ) ;; outer is the source, true if the source is more extreme than the citizen
    and inner * outer > 0          ;; true if citizen and source are on the same side of 0.0
  [ set return-value? true ]
  report return-value?
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; function RxR-> R increasing the first argument by the second argument
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to-report increase [ value change ]
  let ttl -1
  let old-value value
  let diff 0.1 * change / exp ( value * value )
  let result value + diff
  trace ttl ( word "value from " old-value " up to " result " by " diff " derived from " change " " ifelse-value (old-value * result ) < 0 [ "zero crossed" ] [ "" ] )
  if old-value > 0 and result < 0 [ set plus-to-minus plus-to-minus + 1 ]
  if old-value < 0 and result > 0 [ set minus-to-plus minus-to-plus + 1 ]
  report result
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; function RxR-> R decreasing the first argument by the second argument
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to-report decrease [ value change ]
  let ttl -1
  let old-value value
  let diff 0.1 * change / exp ( value * value )
  let result value - diff
  trace ttl ( word "value from " old-value " down to " result " by " diff " derived from " change  " " ifelse-value (old-value * result ) < 0 [ "zero crossed" ] [ "" ] )
  if old-value > 0 and result < 0 [ set plus-to-minus plus-to-minus + 1 ]
  if old-value < 0 and result > 0 [ set minus-to-plus minus-to-plus + 1 ]
  report result
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; function RxRxR -> R changing the first argument in the direction of the
;; second argument by the third argument such that the result is always in R
;; implements B25
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to-report upd-attitude [ changed changer exposure-rank ]
  let ttl -1
  let old-att changed
  let new-att changed
  if changer > 0 [ set new-att increase new-att changer ]
  if changer < 0 [ set new-att decrease new-att abs changer ]
  trace ttl ( word "attitude " old-att " updated to " new-att ", changer: " changer ", rank: " exposure-rank )
  if abs ( old-att - new-att ) > 0.5 [ trace ttl " great change !!! " ]
  report new-att
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; functions to set and get the elements of the memo structure
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to-report sender [ memo ]
  report item 1 memo
end

to-report broadcaster [ memo ]
  report item 3 memo
end

to-report week [ memo ]
  report item 9 memo
end

to-report opinion [ memo ]
  report item 5 memo
end

to-report aspect [ memo ]
  report item 7 memo
end

to-report aspect.num [ memo ]
  let aspects "abcde"
  report position item 7 memo aspects
end

to-report state [ memo ]
  report item 11 memo
end

to-report forgettability [ memo ]
  report item 13 memo
end

to-report set-broadcaster [ memo via ]
  let ttl -1
  trace ttl ( word "memo " memo " will change item 3 to " via )
  set memo replace-item 3 memo via
  report memo
end

to-report set-state [ memo a-state ]
  set memo replace-item 11 memo a-state
  report memo
end

to-report set-forgettability [ memo a-probability ]
  set memo replace-item 13 memo a-probability
  report memo
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; function -> {-1, 0, 1}
;; implements C1-C2
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to-report vote-or-abstain
  let ttl -1
  let current-vote 0
  set probability-to-vote 1 - 1 / exp ( ( interest-in-the-issue + 4 ) * abs ( attitude-position ) )
  let ran random-float 1.0
  if ran < probability-to-vote
  [
    ifelse attitude-position < 0.5
    [ set current-vote -1 ][ set current-vote 1 ]
  ]
  report current-vote
end



;; utilities

to move-to-position ;; moves places and citizens to their positions, leaving a zone on top and at the right free for channels and sources
  setxy ( random ( max-pxcor - min-pxcor - 2 ) + min-pxcor + 2 ) ( random ( max-pycor - min-pycor - 8 ) + min-pycor + 1 )
end

to convert-to-color [ value ]
  let ttl -1
  let old-color color
  let color-limits   ( list -9 -1.2815515655446 -0.841621233572914 -0.524400512708041 -0.2533471031358 0 0.2533471031358 0.524400512708041 0.841621233572914 1.2815515655446 8.20953615160139 )
  let cvalue 0
  ( ifelse
    value < item 1 color-limits [ set cvalue 1 ]
    value < item 2 color-limits [ set cvalue 2 ]
    value < item 3 color-limits [ set cvalue 3 ]
    value < item 4 color-limits [ set cvalue 4 ]
    value < item 5 color-limits [ set cvalue 5 ]
    value < item 6 color-limits [ set cvalue 6 ]
    value < item 7 color-limits [ set cvalue 7 ]
    value < item 8 color-limits [ set cvalue 8 ]
    value < item 9 color-limits [ set cvalue 9 ]
    [ set cvalue 10 ]
  )
  set color approximate-hsb ( cvalue * 30 - 30 ) 100 100
  trace ttl ( word " converted color from " old-color " to " color )
end

to put-color-legend-to-view
  let legend ( list -1.2815515655446 -0.841621233572914 -0.524400512708041 -0.2533471031358 0 0.2533471031358 0.524400512708041 0.841621233572914
                    1.2815515655446 "+inf" )
  let fields [ 1 2 3 4 5 6 7 8 9 10 ]
  ( foreach legend fields
    [ [ leg field ] -> ask patch min-pxcor ( max-pycor - ( field )  ) [ set pcolor approximate-hsb ( field * 30 - 30 ) 100 100 ] ]
  )
  ( foreach legend fields
    [ [ leg field ] -> ask patch ( min-pxcor + 1 ) ( max-pycor - field )
      [
        sprout-dummies 1
        [
          set size 0
          setxy pxcor + 0.7 pycor - 0.7
          set label-color black
          set label ifelse-value ( is-number? leg ) [ precision leg 2 ] [ leg ]
        ]
      ]
    ]
  )
  ask patch ( min-pxcor + 1 ) ( max-pycor )
    [
      sprout-dummies 1
      [
        set size 0
        setxy pxcor + 0.7 pycor - 0.7
        set label-color black
        set label "-inf"
      ]
    ]
end

to setup-temporary-plot-pens
  set-current-plot "positions over time"
  ask citizens
  [
    create-temporary-plot-pen word "citizen " ( who + 1 )
    set-plot-pen-color colour-at-start
  ]
  create-temporary-plot-pen "zero"
  set-plot-pen-color black
end

to do-plot
  set-current-plot "positions over time"
  ask citizens
  [
    set-current-plot-pen word "citizen " ( who + 1 )
    set-plot-pen-color colour-at-start
    plotxy ticks attitude-position
  ]
  set-current-plot-pen "zero"
  set-plot-pen-color black
  plotxy ticks 0.0
end

to document-memos [ memos ]
  file-open logfile-name
  file-show ( word " has the following " length memos " memos in its memory:" )
  foreach memos file-print
  file-print "-------------------------------------------"
end

to trace [ level text ]
  file-open logfile-name
  if level > test-level
  [
    file-show text
    file-flush
  ]
  file-close
end

to trace-citizen [ level text memos memos? ]
  file-open logfile-name
  if level > test-level or witness?
  [
    file-show text
    if memos? [ document-memos memos ]
    file-flush
  ]
  file-close
end

to trace-other [ level text memos memos? ]
  file-open logfile-name
  if level > test-level
  [
    file-show text
    if memos? [ document-memos memos ]
    file-flush
  ]
  file-close
end

;; hide/show button functions
to update-links-shown
  ifelse show-talk? [ show-talk-to ] [ hide-talk-to ]
  ifelse show-exposition? [ show-expo ] [ hide-expo ]
  ifelse show-attractivity? [ show-attr ] [ hide-attr ]
  ifelse show-uses-source? [ show-used ] [ hide-used ]
  ifelse show-goes-to-place? [ show-goes ] [ hide-goes ]
end

to hide-talk-to
  ask talk-to [ hide-link ]
end

to show-talk-to
  ask talk-to [ show-link ]
end

to hide-expo
  ask are-exposed-to [ hide-link ]
end

to show-expo
  ask are-exposed-to [ show-link ]
end

to hide-used
  ask are-used-by [ hide-link ]
end

to show-used
  ask are-used-by [ show-link ]
end

to hide-attr
  ask are-attracted-by [ hide-link ]
end

to show-attr
  ask are-attracted-by [ show-link ]
end

to hide-goes
  ask go-to [ hide-link ]
end

to show-goes
  ask go-to [ show-link ]
end

;; calculates the maximum height of a histogram
to-report rectify-date
  let date-string date-and-time
  let year substring date-string 23 27
  let month substring date-string 19 22
  if month = "Okt" [ set month "Oct" ]
  if month = "Dez" [ set month "Dec" ]
  if month = "Mär" [ set month "Mar" ]
  if month = "Mai" [ set month "May" ]
  let day substring date-string 16 18
  let AMPM substring date-string 13 15
  let hour substring date-string 0 2
  let mmin substring date-string 3 5
  let sec substring date-string 6 12
  let month-number 0
  set month-number ( 1 + ( position month ( list "Jan" "Feb" "Mar" "Apr" "May" "Jun" "Jul" "Aug" "Sep" "Oct" "Nov" "Dec" ) ) )
  ifelse month-number < 10 [ set month ( word "0" month-number ) ] [ set month ( word month-number ) ]
  let num-hour read-from-string hour
  let num-min read-from-string mmin
  let num-sec read-from-string sec
  let time-of-day ( num-sec + 60 * num-min + 3600 * num-hour ) / 3600
  if ( AMPM = "PM" or AMPM = "pm" ) and time-of-day < 12 [ set hour ( word ( num-hour + 12 ) ) ]
  if ( AMPM = "AM" or AMPM = "am" ) and hour = "12" [ set hour "00" ]
  set date-string ( word year "-" month "-" day "-" hour "-" mmin "-" sec )
  report date-string
end

;; converts a number into a string with leading zeroes of the desired number of digits,
;; used for setting the names of the directories for the run results
to-report format [ digits number ]
  let result ( word round number )
  while [ length result < digits ] [ set result ( word "0" result ) ]
  report result
end

to-report odd? [ number ]
  let result true
  if number mod 2 = 1 [ set result false ]
  report result
end

to-report skewness [ data ]
  let skew 0
  if empty? data [ report skew ]
  let this-mean mean data
  let centered-data map [ i -> ( i - this-mean ) ] data
  let sum-sq-x-minus-mean sum map [ i -> i * i ] centered-data
  let std-dev sqrt ( sum-sq-x-minus-mean / length data )
  let zscores map [ i -> i / std-dev ] centered-data
  set skew sum map [ i -> i * i * i ] zscores / length zscores
  report skew
end

to-report excess [ data ]
  let kurtosis 0
  if empty? data [ report kurtosis ]
  let this-mean mean data
  let centered-data map [ i -> ( i - this-mean ) ] data
  let sum-sq-x-minus-mean sum map [ i -> i * i ] centered-data
  let std-dev sqrt ( sum-sq-x-minus-mean / length data )
  let zscores map [ i -> i / std-dev ] centered-data
  set kurtosis sum map [ i -> i * i * i * i ] zscores / length zscores
  report kurtosis - 3
end

;; calculates the maximum height of a histogram
to-report max-freq [ ncat a-list ]
  let n-max-freq 0
  if empty? a-list [ report 0.1 ]
  let a-table map [ value -> ( round ( value * ncat ) ) / ncat ] a-list
  let mode max modes a-table
;;  set n-max-freq reduce
;;    [ [ m n ] -> ifelse-value ( n = mode ) [ m + 1 ] [ m ] ] ( fput 0 a-table )
;;  let result min a-list + ( n-max-freq - 1 ) / ncat * ( max a-list - min a-list )
  report one-of modes a-table
end
@#$#@#$#@
GRAPHICS-WINDOW
810
10
1478
679
-1
-1
20.0
1
10
1
1
1
0
0
0
1
-16
16
-16
16
1
1
1
ticks
30.0

BUTTON
8
10
81
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

SLIDER
18
81
190
114
n-places
n-places
0
100
20.0
5
1
NIL
HORIZONTAL

SLIDER
18
120
190
153
n-citizens
n-citizens
0
500
70.0
10
1
NIL
HORIZONTAL

SLIDER
17
161
189
194
n-channels
n-channels
0
10
3.0
1
1
NIL
HORIZONTAL

SLIDER
17
198
189
231
n-sources
n-sources
0
10
2.0
1
1
NIL
HORIZONTAL

BUTTON
81
10
144
43
NIL
go
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
212
10
412
158
current citizen positions
NIL
NIL
-3.0
3.0
0.0
10.0
true
false
"set-plot-pen-mode 1\nset-histogram-num-bars 20" "if length [ attitude-position ] of citizens > 2 [ set-plot-x-range floor ( 100 * min [ attitude-position ] of citizens ) / 100 ceiling ( 100 * max [ attitude-position ] of citizens ) / 100 ]"
PENS
"default" 1.0 0 -16777216 true "" "histogram [ attitude-position ] of citizens"

BUTTON
144
10
207
43
NIL
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

PLOT
8
483
212
631
memos in memory
NIL
NIL
0.0
25.0
0.0
10.0
true
false
"set-plot-pen-mode 1\nset-histogram-num-bars 25" "if ticks > 0 and count citizens > 0 and max [ length memos-in-memory ] of citizens > 0  [ set-plot-x-range min [ length memos-in-memory ] of citizens max [ length memos-in-memory ] of citizens ]\nif ticks > 0 [ set-plot-y-range 0 5 ]\nauto-plot-on\nif ticks > 0 [ set-plot-x-range 0 max [ length memos-in-memory ] of citizens ]"
PENS
"default" 1.0 0 -16777216 true "" "histogram [ length memos-in-memory ] of citizens"

PLOT
212
483
412
631
predicted vote
NIL
NIL
-2.0
2.0
0.0
10.0
true
false
"set-plot-pen-mode 1\nset-histogram-num-bars 5" "if ticks > 0 [ set-plot-y-range 0 5 ]\nauto-plot-on"
PENS
"default" 1.0 0 -16777216 true "" "if ticks > 0 [ histogram [ vote ] of citizens ]"

PLOT
212
158
412
306
current citizens interest
NIL
NIL
-3.0
3.0
0.0
10.0
true
false
"set-plot-pen-mode 1\nset-histogram-num-bars 20" "if ticks > 0 [ set-plot-y-range 0 5 ]\nauto-plot-on\nif length [ interest-in-the-issue ] of citizens > 2 [ set-plot-x-range floor ( 100 * min [ interest-in-the-issue ] of citizens ) / 100 ceiling ( 100 * max [ interest-in-the-issue ] of citizens ) / 100 ]"
PENS
"default" 1.0 0 -16777216 true "" "histogram [ interest-in-the-issue ] of citizens"

OUTPUT
7
631
679
797
10

MONITOR
1183
742
1349
787
NIL
max-memory-size-used
17
1
11

PLOT
578
11
778
307
positions over time
NIL
NIL
0.0
10.0
-3.0
3.0
true
false
"" ""
PENS

MONITOR
901
693
995
738
NIL
count talk-to
17
1
11

MONITOR
996
693
1154
738
NIL
count are-attracted-by
17
1
11

MONITOR
904
742
1176
787
NIL
count citizens with [ changed-attitude? ]
17
1
11

PLOT
212
306
412
434
attitude changes first
NIL
NIL
-0.5
0.5
0.0
10.0
true
false
"set-plot-pen-mode 1\nset-histogram-num-bars 20" "set-plot-pen-mode 1\nset-histogram-num-bars 20\nif not empty? attitude-changes-first [ set-plot-x-range floor ( 1000 * min attitude-changes-first ) / 1000 ceiling ( 1000 * max attitude-changes-first ) / 1000 ]"
PENS
"default" 1.0 0 -16777216 true "" "if not empty? attitude-changes-first [ histogram attitude-changes-first ]"

PLOT
495
306
695
434
attitude changes second
NIL
NIL
-0.5
0.5
0.0
10.0
true
false
"set-plot-pen-mode 1\nset-histogram-num-bars 20" "set-plot-pen-mode 1\nset-histogram-num-bars 20\nif ticks > 0 [ set-plot-x-range floor ( 1000 * min attitude-changes-second ) / 1000 ceiling ( 1000 * max attitude-changes-second ) / 1000 ]"
PENS
"default" 1.0 0 -16777216 true "" "histogram attitude-changes-second"

MONITOR
212
436
412
481
rng attitude changes first
ifelse-value length attitude-changes-first > 0 [ ( word ( floor ( 1000 * min attitude-changes-first ) / 1000 ) \" : \" ( ceiling ( 1000 * max attitude-changes-first ) / 1000 ) ) ][ \"no changes yet\" ]
17
1
11

MONITOR
495
436
695
481
rng attitude changes second
( word ( floor ( 1000 * min attitude-changes-second ) / 1000 ) \" : \" ( ceiling ( 1000 * max attitude-changes-second ) / 1000 ) )
17
1
11

MONITOR
413
10
493
55
mean
att-mean
4
1
11

MONITOR
413
55
493
100
std dev
standard-deviation [ attitude-position ] of citizens
4
1
11

MONITOR
413
100
493
145
skewness
att-skew
4
1
11

MONITOR
413
158
493
203
mean
mean [ interest-in-the-issue ] of citizens
4
1
11

MONITOR
413
203
493
248
std dev
standard-deviation [ interest-in-the-issue ] of citizens
4
1
11

MONITOR
413
306
493
351
mean
mean attitude-changes-first
4
1
11

MONITOR
413
351
493
396
std dev
standard-deviation attitude-changes-first
4
1
11

MONITOR
698
306
778
351
mean
mean attitude-changes-second
4
1
11

MONITOR
698
351
778
396
std dev
standard-deviation attitude-changes-second
4
1
11

MONITOR
698
396
778
441
skewness
skewness attitude-changes-second
4
1
11

MONITOR
413
396
493
441
skewness
skewness attitude-changes-first
4
1
11

MONITOR
414
249
493
294
skewness
skewness [ interest-in-the-issue ] of citizens
4
1
11

MONITOR
1158
692
1259
737
memos -0+
( list negative-memos neutral-memos positive-memos )
17
1
11

PLOT
414
481
776
631
attitude distribution
NIL
NIL
0.0
10.0
-0.05
0.05
true
true
"" ""
PENS
"mean" 1.0 0 -16777216 true "" "if ticks > 0 [ plot att-mean ]"
"skew" 1.0 0 -2674135 true "" "if ticks > 0 [ plot att-skew ]"
"median" 1.0 0 -13345367 true "" "if ticks > 0 [ plot att-median ]"

SLIDER
18
49
190
82
max-ticks
max-ticks
0
100
100.0
10
1
NIL
HORIZONTAL

MONITOR
495
10
571
55
median
median [ attitude-position ] of citizens
4
1
11

MONITOR
495
55
571
100
mode
max-freq 20 [ attitude-position ] of citizens
4
1
11

MONITOR
495
100
571
145
excess
excess [ attitude-position ] of citizens
4
1
11

SWITCH
17
248
200
281
show-talk?
show-talk?
0
1
-1000

SWITCH
17
281
200
314
show-exposition?
show-exposition?
1
1
-1000

SWITCH
17
314
200
347
show-attractivity?
show-attractivity?
1
1
-1000

SWITCH
17
347
200
380
show-uses-source?
show-uses-source?
1
1
-1000

SWITCH
17
381
200
414
show-goes-to-place?
show-goes-to-place?
1
1
-1000

BUTTON
17
417
200
450
NIL
update-links-shown
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
1264
691
1364
736
+:- / -:+
( list plus-to-minus minus-to-plus )
17
1
11

@#$#@#$#@
## WHAT IS IT?

This model is another attempt to replicate Robert P. Abelson's and Alex Bernstein's "Computer Simulation Model of Community Referendum" (Public Opinion Quarterly vol. XXVII, no. 1, 1965, pp. 93-122. Abelson and Bernstein had "several purposes in mind: to describe the specific features of this particular simulation model, bringing several levels of theory and both experimental and field phenomena to bear upon the total conception; to illustrate the properties of the model by giving some results of a preliminary trial upon artificial, albeit realistic, data; to discuss some of the broad problems that are likely to be encountered in this type of approach; and finally, thus, to elucidate the general character of simulation technique, which seems to offer eventual promise of uniting theories of individual behavior with theories of group behavior." (p. 93)

## HOW IT WORKS

Citizen agents follow a fairly large number of rules. Rules A1 to A22 determine their information processing with respect to information channels, rules B1 to B27 determine how they process information that they receive from other citizens. Other agents represent sources of information (information being packed into "memos") and channels through which these memos are transmitted. Sources and channels do not change their attributes, sources issue new memos in every round, and citizens change their attitudes according the memos they receive, process, accept or reject (but rejected memos remain in their memories unless they are discarded in the early phase of each round when the memos broadcast to all citizens are selected by each citizen according to their attitudes towards sources and channels and according to their interest in the issue as a whole.

Sources can be found in the top row of the view, channels can be found in the third row from above, cititens and the places where they meet to communicate are distributed over the view except its top 10 rows.

The difference between this and an older version (to be found at https://www.comses.net/codebases/a1a0c5b7-f427-4f32-bd70-60c307b61697/releases/1.0.0/) is mainly that attitudes and interest are now normally distributed with mean 0 and variance 1 during setup (instead of uniformly between 0 and 1). This led to further changes in the algorithms to change these two instance variables of the citizens.

The model writes three files. It first creates a directory named `ex-c-s-p-c` (if it does not exist already); the first c is the number of channels, s is the number of sources, p is the number of places, and the second c is the number of citizens. Within this directory a directory is created named `run_YYYY-MM-DD-hh-mm-ss.mss_runn` (mss is the millisecond when the directory is created, runn is the behaviorspace-run-number when this exists, otherwise it is 0000). This directory will contain the files runn-agents.csv and runn-sources.csv (with runn again the run number or 0000) and a log file call log.txt. If the model is run on the GUI, views are stored at the end of every period, named `AbeBer63-nnn.png`, with nnn the number of the period.

## HOW TO USE IT

An ODD of the former version of the model can also be found at https://www.comses.net/codebases/a1a0c5b7-f427-4f32-bd70-60c307b61697/releases/1.0.0/.

An updated version of the ODD will also be published at ComSES as release 2.0.0.

## THINGS TO NOTICE

Links and turtles change their colour according to their most important attributes. These attributes have values between approximately -3 and +3, and the colours are explained in the upper left corner of the view.

## THINGS TO TRY

With the help of the five `show-...?` switches and the `update-links-shown` button, different link types can be made visible or hidden. The button takes effect at the end of each period.

## EXTENDING THE MODEL

It might be interesting to extend the model to more than one issue (which several aspects each) and to replace the strict time schedule with first and second half of a period with a discrete event-oriented scheduling.

## NETLOGO FEATURES

The model makes wide use of NetLogo's link feature.

## RELATED MODELS

The NetLogo Model Library does not have related models. But see below.

## CREDITS AND REFERENCES

The former version can be downloaded from https://www.comses.net/codebases/a1a0c5b7-f427-4f32-bd70-60c307b61697/releases/1.0.0/

The book chapter (also containing an ODD for the former version) is published as Troitzsch, Klaus G. (2021):_ Formal Design Methods and the Relation Between Simulation Models and Theory: A Philosophy of Science Point of View._ In: T. Rudas, G. Péli (eds.), _Pathways Between Social Science and Computational Social Science,_ Computational Social Science Seriess. Cham: Springer Nature Switzerland, S. 21-46 https://doi.org/10.1007/978-3-030-54936-7_2

This version of the model will be presented in a paper to be published in Frontiers in Political Science, Methods and Measurement, Research Topic: _Research Methodologies in Political Science: The Challenge of AI_:Troitzsch, Klaus G. (2023, to appear): _Early forerunners of computational political science_, See https://www.frontiersin.org/research-topics/52428/research-methodologies-in-political-science-the-challenge-of-ai
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

coin heads
false
0
Circle -7500403 true true 15 15 270
Circle -16777216 false false 22 21 256
Line -16777216 false 165 180 192 196
Line -16777216 false 42 140 83 140
Line -16777216 false 37 151 91 151
Line -16777216 false 218 167 265 167
Polygon -16777216 false false 148 265 75 229 86 207 113 191 120 175 109 162 109 136 86 124 137 96 176 93 210 108 222 125 203 157 204 174 190 191 232 230
Polygon -16777216 false false 212 142 182 128 154 132 140 152 149 162 144 182 167 204 187 206 193 193 190 189 202 174 193 158 202 175 204 158
Line -16777216 false 164 154 182 152
Line -16777216 false 193 152 202 153
Polygon -16777216 false false 60 75 75 90 90 75 105 75 90 45 105 45 120 60 135 60 135 45 120 45 105 45 135 30 165 30 195 45 210 60 225 75 240 75 225 75 210 90 225 75 225 60 210 60 195 75 210 60 195 45 180 45 180 60 180 45 165 60 150 60 150 45 165 45 150 45 150 30 135 30 120 60 105 75

computer workstation
false
0
Rectangle -7500403 true true 60 45 240 180
Polygon -7500403 true true 90 180 105 195 135 195 135 210 165 210 165 195 195 195 210 180
Rectangle -16777216 true false 75 60 225 165
Rectangle -7500403 true true 45 210 255 255
Rectangle -10899396 true false 249 223 237 217
Line -16777216 false 60 225 120 225

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
NetLogo 6.3.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="experiment-4-7-50-500" repetitions="10" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>mean [ attitude-position ] of citizens</metric>
    <metric>standard-deviation [ attitude-position ] of citizens</metric>
    <metric>skewness [ attitude-position ] of citizens</metric>
    <enumeratedValueSet variable="n-channels">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-sources">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-places">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-citizens">
      <value value="500"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment-2-5-50-500" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="n-channels">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-sources">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-places">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-citizens">
      <value value="500"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment-7-3-50-500" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="n-channels">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-sources">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-places">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-citizens">
      <value value="500"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment-1-1-50-500" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="n-channels">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-sources">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-places">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-citizens">
      <value value="500"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment-9-9-50-500" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="n-channels">
      <value value="9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-sources">
      <value value="9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-places">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-citizens">
      <value value="500"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment-c-s-50-500" repetitions="2" sequentialRunOrder="false" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="n-channels">
      <value value="1"/>
      <value value="3"/>
      <value value="5"/>
      <value value="7"/>
      <value value="9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-sources">
      <value value="1"/>
      <value value="3"/>
      <value value="5"/>
      <value value="7"/>
      <value value="9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-places">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-citizens">
      <value value="500"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment-c3-s3-50-500" repetitions="6" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>mean [ attitude-position ] of citizens</metric>
    <metric>standard-deviation [ attitude-position ] of citizens</metric>
    <metric>skewness [ attitude-position ] of citizens</metric>
    <metric>positive-memos</metric>
    <metric>neutral-memos</metric>
    <metric>negative-memos</metric>
    <metric>[ attitude-position ] of sources</metric>
    <enumeratedValueSet variable="n-channels">
      <value value="1"/>
      <value value="4"/>
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-sources">
      <value value="1"/>
      <value value="4"/>
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-places">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-citizens">
      <value value="500"/>
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
