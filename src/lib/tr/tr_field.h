#ifndef TR_FIELD_H
#define TR_FIELD_H


double tr_field(double x);

/**
 * @brief Implement the classical FIELD4 algorithm.
 *
 * @param f The function which defines the diffeerential equation.
 * @param y The initial point to start from.
 * @param h The step.
 * @return The approximated solution.
 *
 * This function returns the solution of dz(x)/dx = @p f(z(x)) on
 * [0,@p h], with z(0) = @p y, that is, it returns z(@p h), using the 4th
 * order Runge-Kutta algorithm.
 */
void tr_field_data_set(double Vmax, double K, double kd, double beta, unsigned int forward);


#endif /* TR_FIELD_H */
