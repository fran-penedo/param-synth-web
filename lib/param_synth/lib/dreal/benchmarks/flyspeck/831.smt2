(set-logic QF_NRA)
(declare-fun x1 () Real)
(declare-fun x2 () Real)
(declare-fun x3 () Real)
(declare-fun x4 () Real)
(declare-fun x5 () Real)
(declare-fun x6 () Real)
(assert (<= 5.5696 x1))
(assert (<= x1 6.3504))
(assert (<= 4.0 x2))
(assert (<= x2 4.7524))
(assert (<= 4.0 x3))
(assert (<= x3 4.7524))
(assert (<= 4.0 x4))
(assert (<= x4 5.0625))
(assert (<= 4.0 x5))
(assert (<= x5 6.3504))
(assert (<= 4.0 x6))
(assert (<= x6 6.3504))
(assert (not (< (+ (* 1.0 (- 1.542)) (+ (* (^ x1 0.5) (- 0.490439)) (+ (* 1.0 (* 0.490439 2.36)) (+ (* (^ x2 0.5) 0.318125) (+ (* 1.0 (* 0.318125 (- 2.0))) (+ (* (^ x3 0.5) 0.32468) (+ (* 1.0 (* 0.32468 (- 2.0))) (+ (* (^ x4 0.5) (- 0.740079)) (+ (* 1.0 (* 0.740079 2.25)) (+ (* (^ x5 0.5) 0.178868) (+ (* 1.0 (* 0.178868 (- 2.0))) (+ (* (^ x6 0.5) 0.205819) (+ (* 1.0 (* 0.205819 (- 2.0))) (+ (/ 3.14159265 2.0) (arctan2 (^ (* 4.0 (* x1 (+ (* x1 (* x4 (+ (- x1) (+ x2 (+ (- x3 x4) (+ x5 x6)))))) (+ (* x2 (* x5 (+ (- x1 x2) (+ x3 (+ (- x4 x5) x6))))) (- (- (- (- (* x3 (* x6 (+ x1 (+ (- x2 x3) (+ x4 (- x5 x6)))))) (* x2 (* x3 x4))) (* x1 (* x3 x5))) (* x1 (* x2 x6))) (* x4 (* x5 x6))))))) 0.5) (- (+ (- (* (- x2) x3) (* x1 x4)) (+ (* x2 x5) (+ (- (* x3 x6) (* x5 x6)) (* x1 (+ (- x1) (+ x2 (+ (- x3 x4) (+ x5 x6)))))))))))))))))))))))) 0.0)))
(check-sat)
(exit)
