;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Routines

(in-package #:com.selwynsimsek.binary-linear-algebra)

;;; Useful macros

(defmacro with-return-binary-array (name dimensions &body rest)
  (assert (symbolp name))
  `(let ((,name (make-array (list ,@dimensions) :element-type 'bit :initial-element 0)))
     ,@rest
     ,name))

;;; Creating vectors and matrices

(defun binary-zero-matrix (m n)
  (make-array (list m n)
              :element-type 'bit
              :initial-element 0))

(defun binary-ones-matrix (m n)
  (make-array (list m n)
              :element-type 'bit
              :initial-element 1))

(defun binary-identity-matrix (m &optional (n m))
  (with-return-binary-array eye (m n)
    (loop for i from 0 below (min n m) do
      (setf (aref eye i i) 1))))

(defun binary-scalar-vector (bit)
  (make-array (list 1) :element-type 'bit :initial-element bit))

(defun binary-scalar-matrix (bit)
  (make-array (list 1 1) :element-type 'bit :initial-element bit))

(defun binary-zero-vector (m)
  (make-array (list m) :element-type 'bit :initial-element 0))

(defun binary-ones-vector (m)
  (make-array (list m) :element-type 'bit :initial-element 1))

(defun binary-one-hot-vector (i n)
  (assert (and (<= 0 i) (< i n)))
  (let ((return-vector (binary-zero-vector n)))
    (setf (aref return-vector i) 1)
    return-vector))

(defun binary-pad-vector (vector from-right &optional (from-left 0))
  (concatenate 'bit-vector (b0v from-left) vector (b0v from-right)))

(defun binary-vector-as-row-vector (vector)
  (with-return-binary-array a (1 (length vector))
    (loop for i from 0 below (length vector) do (setf (aref a 0 i) (aref vector i)))))

(defun binary-vector-as-column-vector (vector)
  (with-return-binary-array a ((length vector) 1)
    (loop for i from 0 below (length vector) do (setf (aref a i 0) (aref vector i)))))

;;; Equality
(defun binary-matrix-matrix-equal-p (matrix-1 matrix-2)
  (destructuring-bind (m1y m1x) (array-dimensions matrix-1)
    (destructuring-bind (m2y m2x) (array-dimensions matrix-2)
      (and (= m1y m2y)
           (= m1x m2x)
           (loop for y below m1y
                 always (loop for x below m1x
                              always (= (aref matrix-1 y x)
                                        (aref matrix-2 y x))))))))

(defun binary-vector-vector-equal-p (v1 v2)
  (and (= (length v1) (length v2))
       (loop for i from 0 below (length v1) always (= (aref v1 i) (aref v2 i)))))

(defun mod2 (n) (mod n 2))

;; Symplectic properties

(defun binary-symplectic-matrix-size (matrix)
  (assert (= (array-rank matrix) 2))
  (assert (= (array-dimension matrix 0) (array-dimension matrix 1)))
  (/ (array-dimension matrix 0) 2))

(defun binary-symplectic-vector-size (vector)
  (assert (evenp (length vector)))
  (/ (length vector) 2))

;; Inner products

(defun binary-vector-inner-product (v1 v2)
                                        ; TODO Can improve this
  (assert (= (length v1) (length v2)))
  (mod2 (loop for i from 0 below (length v1) summing (* (aref v1 i) (aref v2 i)))))

(defun binary-vector-symplectic-inner-product (v1 v2)
  (let ((n (binary-symplectic-vector-size v1)))
    (mod2 (loop for i from 0 below n summing (+ (* (aref v1 i) (aref v2 (+ i n)))
                                                (* (aref v1 (+ i n)) (aref v2 i)))))))

(defun popcount (v)
  (declare (type bit-vector v))
  (loop for bit across v count (not (zerop bit))))

;; Outer products

(defun binary-vector-outer-product (v1 v2)
  (let ((return-matrix (binary-zero-matrix (length v1) (length v2))))
    (loop for i from 0 below (length v1) do
          (loop for j from 0 below (length v2) do
                (setf (aref return-matrix i j) (* (aref v1 i) (aref v2 j)))))
    return-matrix))

;;; Matrix operations

(defun binary-transpose (matrix)
  "Returns the bt of a binary matrix."
  (destructuring-bind (m n) (array-dimensions matrix)
    (let ((result (binary-zero-matrix m n)))
      (loop for i from 0 below m
            do (loop for j from 0 below n
                     do (setf (aref result j i) (aref matrix i j))))
      result)))

(defun binary-matrix-diagonal-vector (matrix)
  "Returns the diagonal vector of a matrix."
  (destructuring-bind (m n) (array-dimensions matrix)
    (let ((return-vector (binary-zero-vector (min m n))))
      (loop for i from 0 below (min m n) do
        (setf (aref return-vector i) (aref matrix i i)))
      return-vector)))

(defun vector-iota (n)
  (loop with vector = (make-array (list n) :element-type 'integer :adjustable nil :displaced-to nil)
        for i from 0 below n do (setf (aref vector i) i) finally (return vector)))


(defun two-arg-binary-matrix-matrix-product (matrix-1 matrix-2)
  (bind (((m k) (array-dimensions matrix-1))
         ((k2 n) (array-dimensions matrix-2)))
    (assert (= k k2))
    (with-return-binary-array product (m n)
      (loop for i from 0 below m do
            (loop for j from 0 below n do
                  (setf (aref product i j) (mod2 (loop for s from 0 below k summing
                                                                            (* (aref matrix-1 i s)
                                                                               (aref matrix-2 s j))))))))))

(defun binary-matrix-matrix-product (&rest matrices)
  (reduce #'two-arg-binary-matrix-matrix-product matrices))

(defun binary-matrix-vector-product (matrix vector)
  (bind (((m n) (array-dimensions matrix)))
    (assert (= n (length vector)))
    (with-return-binary-array b (m)
      (loop for i from 0 below m do (setf (aref b i) (mod2
                                                      (loop for j from 0 below n
                                                            summing
                                                            (* (aref matrix i j)
                                                               (aref vector j)))))))))

(defun two-arg-binary-matrix-matrix-sum (matrix-1 matrix-2))

(defun binary-matrix-matrix-sum (&rest matrices))

(defun two-arg-binary-vector-vector-sum (vector-1 vector-2)
  (bit-xor vector-1 vector-2))
(defun binary-vector-vector-sum (&rest vectors)
  (reduce #'two-arg-binary-vector-vector-sum vectors))

(defun invert-binary-matrix (matrix)
  (error "not implemented yet"))

;;; Destructive operations (that alter the arguments)

;;; Random routines

(defun random-bit () (random 2))

(defun random-binary-vector (n)
  (let ((return-vector (binary-zero-vector n)))
    (loop for i from 0 below n do (setf (aref return-vector i) (random-bit)))
    return-vector))

(defun random-binary-matrix (n &optional (m n))
  (let ((return-matrix (binary-zero-matrix n m)))
    (loop for i from 0 below n do
      (loop for j from 0 below m do
        (setf (aref return-matrix i j) (random-bit))))
    return-matrix))

(defun random-binary-diagonal-matrix (n &optional (m n))
  (let ((random-matrix (binary-zero-matrix n m)))
    (loop for i from 0 below (min m n) do (setf (aref random-matrix i i) (random-bit)))
    random-matrix))

(defun random-binary-upper-triangular-unit-matrix (m &optional (n m))
  (let ((return-matrix (binary-identity-matrix m n)))
    (loop for i from 0 below m do
      (loop for j from (1+ i) below n do
        (setf (aref return-matrix i j) (random-bit))))
    return-matrix))

(defun random-binary-lower-triangular-unit-matrix (m &optional (n m))
  (let ((return-matrix (binary-identity-matrix m n)))
    (loop for i from 0 below m do
      (loop for j from 0 below (min i n) do
        (setf (aref return-matrix i j) (random-bit))))
    return-matrix))

(defun random-binary-upper-triangular-matrix (m &optional (n m))
  (let ((return-matrix (binary-zero-matrix m n)))
    (loop for i from 0 below m do
      (loop for j from i below n do
        (setf (aref return-matrix i j) (random-bit))))
    return-matrix))

(defun random-binary-lower-triangular-matrix (m &optional (n m))
  (let ((return-matrix (binary-zero-matrix m n)))
    (loop for i from 0 below m do
      (loop for j from 0 upto i do
        (setf (aref return-matrix i j) (random-bit))))
    return-matrix))
(defun as-permutation-matrix (permutation &optional (sub-one t))
  (let* ((vector (cl-permutation::perm.rep permutation))
         (n (1- (length vector)))
         (matrix (binary-zero-matrix n n)))
    (loop for i from 0 below n do (setf (aref matrix (funcall (if sub-one #'1- #'identity)
                                                              (aref vector (1+ i))) i) 1))
    matrix))

(defun mallows-distribution (n &optional (random-state *random-state*))
  "Samples a permutation in $S_n$ with the Mallows distribution.
   Algorithm 3 in Appendix A of https://doi.org/10.1109/TIT.2021.3081415."
  (let ((a (vector-iota n)) ; step 1
        (s (cl-permutation:perm-identity n)))
    (loop for i from 1 upto n ; step 2
          do (let* ((m (length a)) ; step 3
                    (random (random 1.0d0 random-state))
                    (k (ceiling (log (1+ (* random (1- (expt 2 m)))) 2))) ; step 4
                    (j (aref a (- m k)))) ; step 5
               (setf (aref (cl-permutation::perm.rep s) i) j) ; step 6
               (setf a (remove j a)))) ; step 7
    s)) ; how to test for this?
(defun random-binary-invertible-matrix (n)
  "A binary random invertible matrix. Taken from Algorithm 4 of https://doi.org/10.1109/TIT.2021.3081415."
  (let ((s (mallows-distribution n)) ; step 1
        (l (bim n))
        (r (bim n))) ; step 2
    (loop for j from 1 upto n do ; step 3
                                 (loop for i from (1+ j) upto n do ; step 4
                                                                   (let ((b (random-bit)))
                                                                     (setf (aref l (1- i) (1- j)) b) ; step 5
                                                                     (when (< (aref (cl-permutation::perm.rep s) i)
                                                                              (aref (cl-permutation::perm.rep s) j)) ; step 6
                                                                       (let ((b (random 2)))
                                                                         (setf (aref r (1- i) (1- j)) b)))))) ; step 7
    (bmm* l (bmm* (as-permutation-matrix s nil) r))))

(defun random-binary-symplectic-matrix (n)
  (error "not implemented yet"))

;;;; Decompositions

(defun plfq-decomposition (matrix)
  "Returns the PLFQ decomposition (together with the rank)."
  (error "not implemented yet"))
