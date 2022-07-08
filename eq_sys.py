{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 36,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "System:\n",
      "0.32*x1 + -0.05*x2 + 0.11*x3 + -0.08*x4 = -2.15\n",
      "0.11*x1 + 0.16*x2 + -0.28*x3 + -0.06*x4 = 0.83\n",
      "0.08*x1 + -0.15*x2 + 0.0*x3 + 0.12*x4 = -1.16\n",
      "-0.21*x1 + 0.13*x2 + -0.27*x3 + 0.0*x4 = -0.44\n",
      "\n",
      "Current solution: [0. 0. 0. 0.]\n",
      "Current solution: [-6.71875  5.1875      -inf     -inf]\n",
      "Current solution: [ nan -inf  inf -inf]\n",
      "Current solution: [-inf  nan  nan  nan]\n",
      "Current solution: [nan nan nan nan]\n",
      "Current solution: [nan nan nan nan]\n",
      "Current solution: [nan nan nan nan]\n",
      "Current solution: [nan nan nan nan]\n",
      "Current solution: [nan nan nan nan]\n",
      "Current solution: [nan nan nan nan]\n",
      "Solution:\n",
      "[nan nan nan nan]\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "/home/daria/.local/lib/python3.6/site-packages/ipykernel_launcher.py:33: RuntimeWarning: divide by zero encountered in double_scalars\n"
     ]
    }
   ],
   "source": [
    "from pprint import pprint\n",
    "import numpy as np\n",
    "\n",
    "ITERATION_LIMIT = 10\n",
    "\n",
    "A = np.array([[0.32, -0.05, 0.11, -0.08],\n",
    "            [0.11, 0.16, -0.28, -0.06],\n",
    "            [0.08, -0.15, 0.00, 0.12],\n",
    "            [-0.21, 0.13, -0.27,0.00]\n",
    "            ])\n",
    "b = np.array([-2.15, 0.83, -1.16, -0.44])\n",
    "\n",
    "print(\"System:\")\n",
    "for i in range(A.shape[0]):\n",
    "    row = [\"{}*x{}\".format(A[i, j], j + 1) for j in range(A.shape[1])]\n",
    "    print(\" + \".join(row), \"=\", b[i])\n",
    "print()\n",
    "\n",
    "x = np.zeros_like(b)\n",
    "\n",
    "for it_count in range(ITERATION_LIMIT):\n",
    "    print(\"Current solution:\", x)\n",
    "    x_new = np.zeros_like(x)\n",
    "\n",
    "    for i in range(A.shape[0]):\n",
    "        s1 = np.dot(A[i, :i], x[:i])\n",
    "        s2 = np.dot(A[i, i + 1:], x[i + 1:])\n",
    "        x_new[i] = (b[i] - s1 - s2) / A[i, i]\n",
    "\n",
    "    if np.allclose(x, x_new, atol=1e-10, rtol=0.):\n",
    "        break\n",
    "\n",
    "    x = x_new\n",
    "\n",
    "print(\"Solution:\")\n",
    "print(x)\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.6.7"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 2
}
