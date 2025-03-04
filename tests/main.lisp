(defpackage com.selwynsimsek.binary-linear-algebra/tests
  (:use :cl
        :com.selwynsimsek.binary-linear-algebra
        :fiveam)
  (:shadowing-import-from #:zr-utils #:define-function)
  (:shadowing-import-from #:metabang-bind #:bind))
(in-package :com.selwynsimsek.binary-linear-algebra/tests)

(def-suite binary-linear-algebra
  :description "Test all of binary-linear-algebra")

(in-suite binary-linear-algebra)

(test test (is (= 2 (+ 1 1))))

;;; Generators for random testing

(defun gen-binary-vector (n) (lambda () (random-binary-vector n)))
(defun gen-binary-matrix (n &optional (m n)) (lambda () (random-binary-matrix n m)))
(defun gen-binary-invertible-matrix (n) (lambda () (random-binary-invertible-matrix n)))

(test vector-sum-zero-p
    (for-all ((n (gen-integer :max 20 :min 0)))
      (for-all ((a (gen-binary-vector n)))
        (is (bvv= (b0v n) (bvv+ a a))))))

(test vector-equal-p
  (for-all ((n (gen-integer :max 20 :min 0)))
    (for-all ((b (gen-binary-vector n))
              (c (gen-binary-vector n)))
      (is (eq (bvv= b c) (equalp b c))))))

(test matrix-equal-p
  (for-all ((m (gen-integer :max 20 :min 0))
            (n (gen-integer :max 20 :min 0)))
    (for-all ((b (gen-binary-matrix m n))
              (c (gen-binary-matrix m n)))
      (is (eq (bmm= b c) (equalp b c))))))

(test multiplication-associative-p
  (for-all ((m (gen-integer :max 5 :min 0))
            (n (gen-integer :max 5 :min 0))
            (k (gen-integer :max 5 :min 0))
            (l (gen-integer :max 5 :min 0)))
    (for-all ((a (gen-binary-matrix m n))
              (b (gen-binary-matrix n k))
              (c (gen-binary-matrix k l)))
      (is (bmm= (bmm* a (bmm* b c))
                (bmm* (bmm* a b) c))))))

(test pluq
  (for-all ((n (gen-integer :max 10 :min 0))
            (m (gen-integer :max 10 :min 0)))
    (for-all ((a (gen-binary-matrix m n)))
      (bind (((:values p l f q r) (pluqr a)))
        (is (binary-matrix-permutation-p p))
        (is (binary-matrix-permutation-p q))
        (is (binary-matrix-lower-triangular-unit-p l))
        (is (binary-matrix-upper-triangular-unit-p f))
        (is (bmm= a (bmm* p l f q)))))))

(test inverse
  (for-all ((n (gen-integer :max 10 :min 0)))
    (for-all ((a (gen-binary-invertible-matrix n)))
      (is (bmm= (binary-identity-matrix n)
                (bmm* a (invert-binary-matrix a))))
      (is (bmm= (binary-identity-matrix n)
                (bmm* (invert-binary-matrix a) a))))))


(test test-2 (is (= 4 (+ 2 2))))
