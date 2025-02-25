;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Routines

(in-package #:com.selwynsimsek.binary-linear-algebra)

;;; Types
(deftype binary-matrix () '(simple-array bit (* *)))
(deftype binary-vector () '(simple-array bit (*)))
;;; Useful macros

(defmacro with-return-binary-array (name dimensions &body rest)
  (assert (symbolp name))
  `(let ((,name (make-array (list ,@dimensions) :element-type 'bit :initial-element 0)))
     ,@rest
     ,name))
(defmacro with-block-decomposition (a b c d (matrix upper-block-width
                                             &optional (upper-block-height upper-block-width))
                                    &body body)
  ;`(progn ,upper-block-width ,upper-block-height ,@body)
  (assert (typep a 'symbol))
  (assert (typep b 'symbol))

  (assert (typep c 'symbol))
  (assert (typep d 'symbol))
  (let ((matrix-gensym (gensym "MATRIX"))
        (width-gensym (gensym "WIDTH"))
        (height-gensym (gensym "HEIGHT"))
        (row-1 (gensym "ROW-1"))
        (row-2 (gensym "ROW-2")))
    `(let ((,matrix-gensym ,matrix)
           (,width-gensym ,upper-block-width)
           (,height-gensym ,upper-block-height))
       (multiple-value-bind (,row-1 ,row-2)
           (split-rowwise ,matrix-gensym ,height-gensym)
         (multiple-value-bind (,(if a a (gensym "A")) ,(if b b (gensym "B")))
             (split-columnwise ,row-1 ,width-gensym)
           (multiple-value-bind (,(if c c (gensym "C")) ,(if d d (gensym "D")))
               (split-columnwise ,row-2 ,width-gensym)
             ,@body))))))
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

;; Properties

(defun binary-vector-p (v)
  (and (arrayp v)
       (= 1 (array-rank v))
       (eq 'bit (array-element-type v))))

(defun binary-matrix-p (v)
  (and (arrayp v)
       (= 2 (array-rank v))
       (eq 'bit (array-element-type v))))

(defun binary-matrix-square-p (m)
  (and (binary-matrix-p m)
       (= (array-dimension m 0) (array-dimension m 1))))

(defun binary-matrix-permutation-p (m)
  (and (binary-matrix-square-p m)
       (bind ((n (array-dimension m 0)))
         (loop for i from 0 below n always
                                    (and (= 1 (loop for j from 0 below n counting (not (zerop (aref m i j)))))
                                         (= 1 (loop for j from 0 below n counting (not (zerop (aref m j i))))))))))

(defun binary-matrix-lower-triangular-p (a)
  (and (binary-matrix-p a)
       (bind (((m n) (array-dimensions a)))
         (loop for i from 0 below m always
               (loop for j from (1+ i) below n always (zerop (aref a i j)))))))

(defun binary-matrix-upper-triangular-p (a)
  (and (binary-matrix-p a)
       (bind (((m n) (array-dimensions a)))
         (loop for i from 0 below m always
                                    (loop for j from 0 below (min i n) always (zerop (aref a i j)))))))

(defun binary-matrix-unit-diagonal-p (m)
  (and (binary-matrix-p m)
       (loop for i from 0 below (min (array-dimension m 0) (array-dimension m 1))
             always (= 1 (aref m i i)))))

(defun binary-matrix-lower-triangular-unit-p (m)
  (and (binary-matrix-unit-diagonal-p m)
       (binary-matrix-lower-triangular-p m)))

(defun binary-matrix-upper-triangular-unit-p (m)
  (and (binary-matrix-unit-diagonal-p m)
       (binary-matrix-upper-triangular-p m)))

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

(defun two-arg-binary-matrix-matrix-sum (matrix-1 matrix-2)
  (bind (((m n) (array-dimensions matrix-1))
         ((m1 n1) (array-dimensions matrix-2)))
    (assert (= m m1))
    (assert (= n n1))
    (with-return-binary-array return-array (m n)
      (loop for i from 0 below m do
            (loop for j from 0 below n do
                  (setf (aref return-array i j)
                        (mod2 (+ (aref matrix-1 i j) (aref matrix-2 i j)))))))))

(defun binary-matrix-matrix-sum (&rest matrices)
  (reduce #'two-arg-binary-matrix-matrix-sum matrices))

(defun two-arg-binary-vector-vector-sum (vector-1 vector-2)
  (bit-xor vector-1 vector-2))
(defun binary-vector-vector-sum (&rest vectors)
  (reduce #'two-arg-binary-vector-vector-sum vectors))

(defun trmm (a b)
  "Computes the product of an upper triangular matrix with another one."
  (destructuring-bind (m n) (array-dimensions b)
    (assert (= m (array-dimension a 0) (array-dimension a 1)))
    (if (= m 1)
        (if (zerop (aref a 0 0)) (apply #'binary-zero-matrix (array-dimensions b)) b)
        (if (= m 0)
            (binary-zero-matrix 0 n)
            (let ((split-point (floor (/ m 2))))
              (with-block-decomposition a1 a2 nil a3 (a split-point)
                (multiple-value-bind (b1 b2) (split-rowwise b split-point)
                  (let ((c1 (trmm a1 b1)))
                    (setf c1 (bmm+ c1 (bmm* a2 b2)))
                    (let ((c2 (trmm a3 b2)))
                      (stack-matrices c1 c2))))))))))




(defun trtri (a)
  "Upper triangular matrix inversion"
  (let ((n (array-dimension a 0)))
    (if (<= n 1)
        (bim n)
        (let ((pivot (floor (/ n 2))))
          (with-block-decomposition a1 a2 nil a3 (a pivot)
            (let* ((c1 (trtri a1))
                   (c3 (trtri a3))
                   (c2 (bmm* a2 c3)))
              (setf c2 (trmm c1 c2))
              (block-matrix c1 c2 nil c3)))))))

(defun ltrtri (a)
  "Lower triangular matrix inversion" ; TODO Improve this implementation
  (bt (trtri (bt a))))

(defun invert-binary-matrix (matrix)
  (bind (((:values p l f q) (pluq-decomposition matrix)))
    (bmm* (bt q) (trtri f) (ltrtri l) (bt p))))

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

(defun block-matrix  (a1 a2 a3 a4)
  "Creates a block matrix. $A_2$ or $A_3$ may be nil, in which case zeroes are filled in."
  (destructuring-bind (m1 n1) (array-dimensions a1)
    (destructuring-bind (m4 n4) (array-dimensions a4)
      (let ((return-matrix (binary-zero-matrix (+ m1 m4) (+ n1 n4))))
        (loop for i from 0 below m1 do
          (loop for j from 0 below n1 do
            (setf (aref return-matrix i j) (aref a1 i j))))
        (loop for i from 0 below m4 do
          (loop for j from 0 below n4 do
            (setf (aref return-matrix (+ i m1) (+ j n1)) (aref a4 i j))))
        (when a2
          (loop for i from 0 below m1 do
            (loop for j from 0 below n4 do
              (setf (aref return-matrix i (+ j n1)) (aref a2 i j)))))
        (when a3
          (loop for i from 0 below m4 do
            (loop for j from 0 below n1 do
              (setf (aref return-matrix (+ i m1) j) (aref a3 i j)))))
        return-matrix))))

(defun split-columnwise
    (a &optional (split-index (floor (/ (array-dimension a 1) 2))))
  "Returns two matrices $A_!$ and $A_2$ that are the left and right halves of $A$ split columnwise."
  (destructuring-bind (m n) (array-dimensions a)
    (let* ((a1 (binary-zero-matrix m split-index))
           (a2 (binary-zero-matrix m (- n split-index))))
      (loop for i from 0 below m do
            (loop for j from 0 below split-index do
              (setf (aref a1 i j) (aref a i j))))
      (loop for i from 0 below m do
        (loop for j from split-index below n do
          (setf (aref a2 i (- j split-index)) (aref a i j))))
      (values a1 a2))))

(defun split-rowwise (a r1)
  "Returns two matrices $A_3$ and $A_4$ that are the upper and lower halves of $A$ split rowwise,
   and that $A_3$ has $r_1$ rows."
  (destructuring-bind (m n) (array-dimensions a)
    (let ((a3 (binary-zero-matrix r1 n))
          (a4 (binary-zero-matrix (- m r1) n)))
      (loop for i from 0 below r1 do
        (loop for j from 0 below n do
          (setf (aref a3 i j) (aref a i j))))
      (loop for i from r1 below m do
        (loop for j from 0 below n do
          (setf (aref a4 (- i r1) j) (aref a i j))))
      (values a3 a4))))

(defun ple-decomposition (a)
  "Input: $$A \in \mathbb F_2^{m \times n}$$ , Output: $(P,L,E)$ a PLE decomposition of $$A$$.
   Given any matrix $A \in F_2^{ m \times n}$  of rank $$r$$, there is a PLE decomposition $A = P LE$ where $P$ is a
   $n \times n$ permutation matrix, $L$ is a $m \times r$ lower triangular matrix and $E$ is a $r \times n$ matrix in
   row-echelon form, with unit leading coefficients.
   This method returns $(P,L,E)$ given $A$.
   Implementation taken from arXiv:1204.3735v1."
  (destructuring-bind (m n) (array-dimensions a)
    (if (= n 0)
        (values (bim m) (binary-zero-matrix m 0) (binary-zero-matrix 0 0))
        (if (= n 1)                     ; step 1
            (if (loop for i from 0 below m always (zerop (aref a i 0)))
                (return-from ple-decomposition
                  (values (bim m) (binary-zero-matrix m 0) (binary-zero-matrix 0 1))) ; step 2
                (let* ((j (loop for k from 0 below m when (not (zerop (aref a k 0))) return k))
                       (p (transposition-permutation-matrix 0 j m))) ; step 3
                  (values p (bmm* p a) (bsm 1))))                    ; step 4
            (multiple-value-bind (a1 a2) (split-columnwise a)        ; step 5
              (multiple-value-bind (p1 l1 e1) (ple-decomposition a1) ; step 6
                (let ((r1 (array-dimension e1 0)))
                  (setf a2 (bmm* (bt p1) a2)) ; step 7
                  (multiple-value-bind (a3 a4) (split-rowwise a2 r1)
                    (multiple-value-bind (l1-1 l1-2) (split-rowwise l1 r1)
                      (setf a3 (ltrsm l1-1 a3))          ; step 8
                      (setf a4 (bmm+ a4 (bmm* l1-2 a3))) ; step 9
                      (multiple-value-bind (p2 l2 e2) (ple-decomposition a4) ; step 10
                        (values         ; step 10
                         (bmm* p1 (block-matrix (bim r1) nil nil p2))
                         (block-matrix l1-1 nil (bmm* (bt p2) l1-2) l2)
                         (block-matrix e1 a3 nil e2))))))))))))

(defun upper-trapezoidal-form (e)
  "Returns $(P,F)$ such that $FP=E$ where $P$ is a permutation matrix, $F$ is upper trapezoidal
   and the input $E$ is in upper trapezoidal form."
  (destructuring-bind (m n) (array-dimensions e)
    (let ((permutation-matrix (binary-zero-matrix n n)))
      (loop for i from 0 below m do
        (setf (aref permutation-matrix i (loop for j from 0 below n
                                               when (not (zerop (aref e i j))) return j))
              1))
      (loop for i from m below n do
        (setf (aref permutation-matrix i
                    (loop for j from 0 below n
                          when (loop for k from 0 below n
                                     always (zerop (aref permutation-matrix k j)))
                            return j))
              1))
      (let ((f (bmm* e (bt permutation-matrix))))
        (values permutation-matrix f)))))

(defun pluq-decomposition (a)
  (bind (((:values p l e) (ple-decomposition a))
         ((:values q f) (upper-trapezoidal-form e)))
    (values p l f q (array-dimension e 0))))


(defun transposition-permutation-matrix (i j m)
  "Returns a $m \times m$ permutation matrix interchanging $i$ and $j$."
  (let ((return-matrix (binary-zero-matrix m m)))
    (loop for k from 0 below m do
      (unless (or (= k i) (= k j))
        (setf (aref return-matrix k k) 1)))
    (setf (aref return-matrix i j) 1
          (aref return-matrix j i) 1)
    return-matrix))


(defun ltrsm (a b)
  "Finds $$X$$ such that $$AX=B$$. $A \in \mathbb F_2^{m \times m}$ is non-singular lower triangular,
   $$B \in \mathbb F_2^{m \times n}$$."
  (destructuring-bind (m n) (array-dimensions b)
    (assert (= m (array-dimension a 0) (array-dimension a 1)))
    (if (or (= n 0) (= m 0))
        (binary-zero-matrix m n)
          (if (= m 1)
              (return-from ltrsm b)
              (let ((split-point (floor (/ m 2))))
                (multiple-value-bind (b1 b2) (split-rowwise b split-point)
                  (multiple-value-bind (a-row1 a-row2) (split-rowwise a split-point)
                    (multiple-value-bind (a1 empty-block) (split-columnwise a-row1 split-point)
                      (multiple-value-bind (a2 a3) (split-columnwise a-row2 split-point)
                        (let* ((x1 (ltrsm a1 b1))
                               (x2 (ltrsm a3 (bmm+ b2 (bmm* a2 x1)))))
                          (stack-matrices x1 x2)))))))))))


(defun stack-matrices (x1 x2)
  "Stacks the matrices on top of each other."
  (destructuring-bind (m1 n1) (array-dimensions x1)
    (destructuring-bind (m2 n2) (array-dimensions x2)
      (assert (= n1 n2))
      (let ((result (binary-zero-matrix (+ m1 m2) n1)))
        (loop for i from 0 below m1 do
          (loop for j from 0 below n1 do
            (setf (aref result i j) (aref x1 i j))))
        (loop for i from 0 below m2 do
          (loop for j from 0 below n2 do
            (setf (aref result (+ i m1) j) (aref x2 i j))))
        result))))
