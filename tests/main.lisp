(defpackage binary-linear-algebra/tests/main
  (:use :cl
        :binary-linear-algebra
        :rove))
(in-package :binary-linear-algebra/tests/main)

;; NOTE: To run this test file, execute `(asdf:test-system :binary-linear-algebra)' in your Lisp.

(deftest test-target-1
  (testing "should (= 1 1) to be true"
    (ok (= 1 1))))
