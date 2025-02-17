;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;; Aliases for functions.

(in-package #:com.selwynsimsek.binary-linear-algebra)

(defmacro function-alias (alias function-name)
  (assert (symbolp alias))
  (assert (symbolp function-name))
  `(progn (setf (fdefinition ',alias) (function ,function-name))
          (export '(,alias ,function-name))))

(function-alias bt binary-transpose)
