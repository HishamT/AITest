class Car:
    def __init__(self, brand, model, year):
        self.brand = brand
        self.model = model
        self.year = year
        self.speed = 0

    def accelerate(acceleration):
        self.speed += acceleration

    def brake(self, deceleration):
        if self.speed > deceleration:
            self.speed -= deceleration
        else:
            self.speed = 0

    def get_speed(self):
        return self.speed

car = Car("Nissan", "Maxima", 2017)
print(car.accelerate(200))
