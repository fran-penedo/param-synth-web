(set-logic QF_NRA)
(declare-fun x () Real)
(assert (< 0.0 x))
(assert (< x 50.0))
(assert (= (^ x 0.5) 7))
(check-sat)
(exit)
