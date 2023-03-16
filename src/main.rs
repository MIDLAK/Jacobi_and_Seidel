use std::io::{BufRead, BufReader};
use std::io;
use rand::Rng;
use std::fs::File;

const EPS: f64 = 0.01; //желаемая точность
const MAX_ITER: i32 = 10_000;
const ERROR_EPS: f64 = 10_000.0; //критичное значение ошибки

///Вывод матрицы любой размерности на экран
fn print_matrix(matrix: Vec<Vec<f64>>) {
    for row in &matrix {
        println!("{row:?}");
    }
}

///Наполнение матрицы случайными значениями от -10.00 до 10.99
fn generate_matrix(matrix: &mut Vec<Vec<f64>>) {
    for row in matrix {
        for el in row {
            *el = rand::thread_rng().gen_range(-10.0..=10.0); 
            *el = (*el * 100.0).round()/100.0;
        }
    }
}

///Вычисление корней СЛАУ методом Якоби
///
///# Параметры
///* `matrix` - Mатрица коэффициентов. Последним столбцом распологаются свободные члены;
///* `x` - Начальное приближение.
///
///# Возвращаемые значения
///Функция возвращает кортеж из значений успеха (true или false) и кол-ва итераций.
fn jacobi(matrix: &mut Vec<Vec<f64>>, x: &mut Vec<f64>) -> (bool, i32) {
    let n = matrix.len();
    let mut iters = 0; //подсчёт итераций

    loop {
        iters += 1;

        let x_prev = x.clone();

        let mut max_eps = 0.0;

        //вычисление очередной итерации
        for i in 0..n {

            if matrix[i][i] == 0.0 {
                println!("Ноль на главной диагонали!");
                return (false, iters);
            }

            let coefficient = 1.0 / matrix[i][i];

            //вычисление скобки
            let mut difference = matrix[i][n];
            for j in 0..n {
                if i != j {
                    difference -= matrix[i][j] * x_prev[j];
                }
            }

            //значение x на очередной итерации
            x[i] = coefficient * difference;

            //значение ошибки
            let eps = (x[i] - x_prev[i]).abs();

            //поиск максимальной ошибки на текущей итерации
            if eps >= max_eps {
                max_eps = eps;
            }
        }

        if max_eps <= EPS {
            return (true, iters);
        }

        if iters >= MAX_ITER || max_eps >= ERROR_EPS {
            return (false, iters);
        }
    }
}

///Вычисление корней СЛАУ методом Зейделя
///
///# Параметры
///* `matrix` - Mатрица коэффициентов. Последним столбцом распологаются свободные члены;
///* `x` - Начальное приближение.
///
///# Возвращаемые значения
///Функция возвращает кортеж из значений успеха (true или false) и кол-ва итераций.
fn seidel(matrix: &mut Vec<Vec<f64>>, x: &mut Vec<f64>) -> (bool, i32) {
    let n = matrix.len();
    let mut iters = 0; //подсчёт итераций

    loop {
        iters += 1;

        let x_prev = x.clone();

        let mut max_eps = 0.0;

        //вычисление очередной итерации
        for i in 0..n {

            if matrix[i][i] == 0.0 {
                println!("Ноль на главной диагонали!");
                return (false, iters);
            }

            let coefficient = 1.0 / matrix[i][i];

            //вычисление скобки
            let mut difference = matrix[i][n];
            for j in 0..n {
                if i != j {
                    difference -= matrix[i][j] * x[j];
                }
            }

            //значение x на очередной итерации
            x[i] = coefficient * difference;

            //значение ошибки
            let eps = (x[i] - x_prev[i]).abs();

            //поиск максимальной ошибки на текущей итерации
            if eps >= max_eps {
                max_eps = eps;
            }
        }

        if max_eps <= EPS {
            return (true, iters);
        }

        if iters >= MAX_ITER || max_eps >= ERROR_EPS {
            return (false, iters);
        }
    }
}

/// Чтение матрицы из файла c именем file_name.
/// Файл должен быть расположен в той же дирректории, где и *.rc файл
fn read_matrix_from_file(file_name: &str) -> Vec<Vec<f64>> {
    let f = BufReader::new(File::open(file_name).unwrap());

    let arr: Vec<Vec<f64>> = f.lines()
        .map(|l| l.unwrap().split(char::is_whitespace)
             .map(|number| number.parse().unwrap())
             .collect())
        .collect();

    return arr;
}

fn main() {
    println!("Точность {EPS}");

    let mut count_gen: usize = 0;
    let mut dimension: usize = 0;

    let mut file_name = String::new();
    let mut file_flag = false;
    
    //цикл для получения верного ввода
    loop {
        println!("Введите имя файла (или оставьте поле пустым для генерации): ");

        io::stdin()
            .read_line(&mut file_name)
            .expect("Ошибка ввода!");

        println!("Введите размерность матрицы: ");

        let mut dimension_str = String::new();

        io::stdin()
            .read_line(&mut dimension_str)
            .expect("Ошибка ввода!");

        dimension = match dimension_str.trim().parse() {
            Ok(num) => num,
            Err(_) => continue,
        };

        count_gen = 1;
        
        if file_name.trim().is_empty() {
            println!("Количество генераций: ");

            let mut count_gen_str = String::new();

            io::stdin()
                .read_line(&mut count_gen_str)
                .expect("Ошибка ввода!");
            
            count_gen = match count_gen_str.trim().parse() {
                Ok(num) => num,
                Err(_) => continue,
            };
        } else {
            file_flag = true;
        }

        break;
    }

    let mut iters = 0;
        
    let mut success = 0.0;
    let mut average_iters = 0.0;

    let mut success_seidel = 0.0;
    let mut average_iters_seidel = 0.0;

    loop {
        let mut matrix = vec![vec![0.0f64; dimension+1]; dimension];

        if file_flag {
            let arr = read_matrix_from_file(file_name.trim());
            matrix = arr.clone();
        } else {
            generate_matrix(&mut matrix);
        }
        
        let mut x = vec![0.0f64; dimension];

        let result_seidel = seidel(&mut matrix, &mut x);
        if result_seidel.0 {
            success_seidel += 1.0;
            average_iters_seidel += result_seidel.1 as f64;
        }
        
        x = vec![0.0f64; dimension];

        let result_jacobi = jacobi(&mut matrix, &mut x);
        if result_jacobi.0 {
            success += 1.0;
            average_iters += result_jacobi.1 as f64;
        } 

            
        iters += 1;

        if iters >= count_gen { 
            println!("ЯКОБИ");
            let ratio = 100.0 * success / count_gen as f64;
            let average_iters = average_iters / success as f64;
            println!("Сошлись {success} уравнений из {count_gen} ({ratio:.4} %).");
            println!("Среднее число итераций составило {average_iters:.4}.");

            println!("ЗЕЙДЕЛЬ");
            let ratio_seidel = 100.0 * success_seidel / count_gen as f64;
            let average_iters_seidel = average_iters_seidel / success_seidel as f64;
            println!("Сошлись {success_seidel} уравнений из {count_gen} ({ratio_seidel:.4} %).");
            println!("Среднее число итераций составило {average_iters_seidel:.4}.");

            break;
        }
    }
}
