(defpackage com.selwynsimsek.binary-linear-algebra/tests
  (:use :cl
        :com.selwynsimsek.binary-linear-algebra
        :fiveam)
  (:shadowing-import-from #:zr-utils #:define-function)
  (:shadowing-import-from #:metabang-bind #:bind))
(in-package :com.selwynsimsek.binary-linear-algebra/tests)


#+continuous-integration
(setf fiveam:*debug-on-error* t fiveam:*debug-on-failure* t)

(def-suite binary-linear-algebra
  :description "Test all of binary-linear-algebra")

(in-suite binary-linear-algebra)

;;; Generators for random testing

(defun gen-binary-vector (n) (lambda () (random-binary-vector n)))
(defun gen-binary-matrix (n &optional (m n)) (lambda () (random-binary-matrix n m)))
(defun gen-binary-invertible-matrix (n) (lambda () (random-binary-invertible-matrix n)))
(defun gen-binary-invertible-upper-triangular-unit-matrix (n)
  (lambda () (random-binary-upper-triangular-unit-matrix n)))

(test vector-sum-zero-p
    (for-all* ((n (gen-integer :max 20 :min 0))
               (a (gen-binary-vector n)))
      (is (bvv= (b0v n) (bvv+ a a)))))

(test vector-equal-p
  (for-all* ((n (gen-integer :max 20 :min 0))
             (b (gen-binary-vector n))
             (c (gen-binary-vector n)))
    (is (eq (bvv= b c) (equalp b c)))))

(test matrix-equal-p
  (for-all* ((m (gen-integer :max 20 :min 0))
             (n (gen-integer :max 20 :min 0))
             (b (gen-binary-matrix m n))
             (c (gen-binary-matrix m n)))
    (is (eq (bmm= b c) (equalp b c)))))

(test multiplication-associative-p
  (for-all* ((m (gen-integer :max 5 :min 0))
             (n (gen-integer :max 5 :min 0))
             (k (gen-integer :max 5 :min 0))
             (l (gen-integer :max 5 :min 0))
             (a (gen-binary-matrix m n))
             (b (gen-binary-matrix n k))
             (c (gen-binary-matrix k l)))
    (is (bmm= (bmm* a (bmm* b c))
              (bmm* (bmm* a b) c)))))

(test multiplication-distributive-p
 (for-all* ((m (gen-integer :max 10 :min 0))
            (n (gen-integer :max 10 :min 0))
            (k (gen-integer :max 10 :min 0))
            (a (gen-binary-matrix m n))
            (b (gen-binary-matrix n k))
            (c (gen-binary-matrix n k)))
   (is (bmm= (bmm+ (bmm* a b) (bmm* a c))
             (bmm* a (bmm+ b c))))))

(test pluq
  (for-all* ((n (gen-integer :max 10 :min 0))
             (m (gen-integer :max 10 :min 0))
             (a (gen-binary-matrix m n)))
    (bind (((:values p l f q r) (pluqr a)))
      (is (binary-matrix-permutation-p p))
      (is (binary-matrix-permutation-p q))
      (is (binary-matrix-lower-triangular-unit-p l))
      (is (binary-matrix-upper-triangular-unit-p f))
      (is (bmm= a (bmm* p l f q))))))

(test inverse
  (for-all* ((n (gen-integer :max 10 :min 0))
             (a (gen-binary-invertible-matrix n)))
    (is (bmm= (binary-identity-matrix n)
              (bmm* a (invert-binary-matrix a))))
    (is (bmm= (binary-identity-matrix n)
              (bmm* (invert-binary-matrix a) a)))))

(test invertible-p
      (for-all* ((n (gen-integer :max 10 :min 0))
                 (a (gen-binary-matrix n)))
               (is (if (ignore-errors (ibm a))
                       (binary-matrix-invertible-p a)
                       (not (binary-matrix-invertible-p a))))))

(test trsm
      (for-all* ((n (gen-integer :max 10 :min 0))
                 (m (gen-integer :max 10 :min 0))
                 (a (gen-binary-invertible-upper-triangular-unit-matrix m))
                 (b (gen-binary-matrix m n)))
        (is (bmm= b (bmm* a (trsm a b))))))

(test solve-matrix-system
  (for-all* ((n (gen-integer :max 10 :min 0))
             (m (gen-integer :max 10 :min 0))
             (k (gen-integer :max 10 :min 0))
             (a (gen-binary-matrix n m))
             (b (gen-binary-matrix n k)))
    (let ((x (solve-matrix-system a b)))
      (is (if x
              (bmm= b (bmm* a x))
              t)))))

(test row-echelon-form
  (for-all* ((n (gen-integer :max 10 :min 0))
             (m (gen-integer :max 10 :min 0))
             (a (gen-binary-matrix n m)))
    (bind (((:values e x) (row-echelon-form a)))
      (is (bmm= e (bmm* x a)))
      (is (binary-matrix-invertible-p x))
      (is (row-echelon-p e)))))

(test reduced-row-echelon-form
  (for-all* ((n (gen-integer :max 10 :min 0))
             (m (gen-integer :max 10 :min 0))
             (a (gen-binary-matrix n m)))
    (bind (((:values r y) (reduced-row-echelon-form a)))
      (is (bmm= r (bmm* y a)))
      (is (binary-matrix-invertible-p y))
      (is (reduced-row-echelon-p r)))))

(test binary-symplectic-matrix
  (for-all* ((n (gen-integer :max 10 :min 0)))
    (let ((m (random-binary-symplectic-matrix n)))
      (is (binary-matrix-symplectic-p m)))))
