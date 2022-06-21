# Библиотека схемы разделения секрета на основе матриц Адамара

Схема, реализованная по статье <a href="https://ieeexplore.ieee.org/document/8359980">A Secret Sharing Scheme from Hadamard Matrix</a> с некоторыми изменениями и дополнениями.

### Отличие от статьи
Предполагается, что доли, полученные при разделении секрета, содержат в себе только номер строки матрицы и значения, полученные по данной строке. Хранение матрицы предполагается делегировать дилеру и не раскрывать её вид участникам. Данное изменении используется для достижения свойства совершенности.

### Дополнения
Возможность валидации набора долей или обнаружения злоумышленника.

## Скачивание
```bash
git clone git@github.com:DimartX/hadamard-secret-sharing.git
```

## Получение документации
Библиотека содержит небольшую документацию. Получить её можно следующим образом.

```bash
cargo doc --open
```

## Сборка и тестирование
Сборка библиотеки.

```bash
cargo build
```

Прогонка юнит-тестов.
```bash
cargo test --lib
```

## Использование
Для использования библиотеки в своём проекте необходимо указать в `Cargo.toml` зависимость.

```toml
hadamard_sss = { git = "https://github.com/DimartX/hadamard-secret-sharing.git", version = "0.1.1" }
```

Далее приводится пример из <a href="https://github.com/DimartX/hadamard_sss_example">hadamard_sss_example</a>.

### Пояснение примера `simple_example`
Добыча матрицы отдаётся на откуп пользователю библиотеки. Предполагается хранение матрицы у дилера.

Считываем матрицу Адамара (взята с http://neilsloane.com/hadamard/) из файла.

```rust
let mtx = read_matrix("../matrices/had8.txt").expect("can't read matrix");
```

Создаём экземляр структуры схемы разделения секрета.

```rust
let scheme = HadamardSSS::from(&mtx).expect("can't create scheme");
```

Создаём секрет типа u32 и разделяем его на доли, количество которых на одну меньше размерности используемой матрицы.

```rust
let secret: u32 = 314159265;
let parts = scheme.share(secret).expect("can't share that secret");
```

Демонстрируем результат при попытке восстановить секрет по числу долей, меньшему порогового значения.

```rust
let parts_rec = vec![parts[1], parts[3], parts[6]];
let res_secret = scheme.reconstruct(parts_rec);
match res_secret {
    Ok(res) => println!("Result is {}", res),
    Err(err_msg) => println!("can't reconstruct secret: {}", err_msg),
}
```

Также демонстрируем успешное восстановление секрета.

```rust
let parts_rec = vec![parts[1], parts[3], parts[4], parts[5], parts[6]];
let res_secret = scheme.reconstruct(parts_rec);
match res_secret {
    Ok(res) => println!("Result is {}", res),
    Err(err_msg) => println!("can't reconstruct secret: {}", err_msg),
}
```

Результат выполнения примера `simple_example`:

```
Matrix was read:
[[1, 1, 1, 1, 1, 1, 1, 1],
 [1, -1, 1, -1, 1, -1, 1, -1],
 [1, 1, -1, -1, 1, 1, -1, -1],
 [1, -1, -1, 1, 1, -1, -1, 1],
 [1, 1, 1, 1, -1, -1, -1, -1],
 [1, -1, 1, -1, -1, 1, -1, 1],
 [1, 1, -1, -1, -1, -1, 1, 1],
 [1, -1, -1, 1, -1, 1, 1, -1]]
The initial secret is 314159265
Parts:
(0, 1480438321)
(1, 446694049)
(2, 114636929)
(3, 415081689)
(4, 1203303429)
(5, 2973513135)
(6, 2746864041)
Trying to reconstruct by 3 random values:
3 is less than threshold 5 parties
can't reconstruct secret: less than threshold parties
Trying to reconstruct by 5 random values:
Result is 314159265
```
