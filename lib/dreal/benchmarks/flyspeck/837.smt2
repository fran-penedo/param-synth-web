(set-logic QF_NRA)
(declare-fun x1 () Real)
(declare-fun x2 () Real)
(declare-fun x3 () Real)
(declare-fun x4 () Real)
(declare-fun x5 () Real)
(declare-fun x6 () Real)
(assert (<= 4.0 x1))
(assert (<= x1 6.3504))
(assert (<= 4.0 x2))
(assert (<= x2 6.3504))
(assert (<= 4.0 x3))
(assert (<= x3 6.3504))
(assert (<= 6.3504 x4))
(assert (<= x4 8.0))
(assert (<= 4.0 x5))
(assert (<= x5 6.3504))
(assert (<= 4.0 x6))
(assert (<= x6 6.3504))
(assert (not (< (+ (* 1.0 (- 1.08346)) (+ (* 1.0 (* 0.288794 (- 2.0))) (+ (* (^ x1 0.5) 0.288794) (+ (* 1.0 (* 0.292829 2.0)) (+ (* (^ x2 0.5) (- 0.292829)) (+ (* 1.0 (* 0.036457 (- 2.0))) (+ (* (^ x3 0.5) 0.036457) (+ (* 1.0 (* 0.348796 (- 2.52))) (+ (* (^ x4 0.5) 0.348796) (+ (* 1.0 (* 0.762602 2.0)) (+ (* (^ x5 0.5) (- 0.762602)) (+ (* 1.0 (* 0.112679 2.0)) (+ (* (^ x6 0.5) (- 0.112679)) (+ (+ (/ 3.14159265 2.0) (arctan2 (^ (* 4.0 (* x2 (+ (* x2 (* x5 (+ (- x2) (+ x1 (+ (- x3 x5) (+ x4 x6)))))) (+ (* x1 (* x4 (+ (- x2 x1) (+ x3 (+ (- x5 x4) x6))))) (- (- (- (- (* x3 (* x6 (+ x2 (+ (- x1 x3) (+ x5 (- x4 x6)))))) (* x1 (* x3 x5))) (* x2 (* x3 x4))) (* x2 (* x1 x6))) (* x5 (* x4 x6))))))) 0.5) (- (+ (- (* (- x1) x3) (* x2 x5)) (+ (* x1 x4) (+ (- (* x3 x6) (* x4 x6)) (* x2 (+ (- x2) (+ x1 (+ (- x3 x5) (+ x4 x6))))))))))) (* 1.0 (- 0.0031)))))))))))))))) 0.0)))
(check-sat)
(exit)
