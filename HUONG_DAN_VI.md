# Hướng dẫn ngắn

## Chạy bài

Giải nén ZIP, mở terminal trong thư mục chứa mã nguồn rồi chạy:

```bash
python3 -m pip install -r requirements.txt
python3 ellipse_distance.py --seed 42
```

Bấm chuột vào cửa sổ đồ thị, sau đó nhấn **Enter** hoặc **Space** để tìm hai điểm.
**N** tạo ví dụ mới, **S** lưu ảnh, **Q** đóng cửa sổ.
Trên Windows có thể thay `python3` bằng `python`.

## Hiểu thuật toán để giải thích với thầy

1. Một ellipse được xác định bởi tâm, hai bán trục và góc xoay. Chương trình
   tạo ngẫu nhiên các giá trị này cùng một điểm truy vấn (u0,v0).
2. Ta cần tìm điểm nằm trên ĐƯỜNG BIÊN ellipse có khoảng cách nhỏ nhất/lớn
   nhất tới điểm truy vấn. Nếu điểm truy vấn ở trong ellipse, khoảng cách
   nhỏ nhất tới đường biên thường vẫn lớn hơn 0.
3. Tối ưu bình phương khoảng cách cho cùng kết quả với tối ưu khoảng cách,
   nhưng đạo hàm đơn giản hơn vì không có căn bậc hai.
4. Dùng nhân tử Lagrange để kết hợp hàm khoảng cách và ràng buộc ellipse.
   Cho gradient bằng 0 thu được một hệ tuyến tính phụ thuộc lambda.
5. Thế nghiệm của hệ vào phương trình ellipse thu được đa thức bậc bốn.
   Tìm các nghiệm thực, dựng lại các điểm ứng viên rồi so sánh khoảng cách.
6. Mã nguồn chuyển về hệ tọa độ có tâm ellipse ở gốc và trục ellipse trùng
   trục tọa độ trước khi áp dụng công thức. Đây vẫn là phương pháp trong PDF.
7. `numpy.roots` giải đa thức bằng số thực dấu phẩy động; chương trình không
   dùng thuật toán tối ưu số như gradient descent hay `scipy.optimize`.
8. Điểm trên trục, điểm ở tâm và đường tròn được xử lý riêng vì phép nghịch
   đảo ma trận trong tài liệu có thể không tồn tại ở nghiệm đúng.

Xem README.md để có đầy đủ công thức và liên hệ giữa công thức PDF với code.

## Đối chiếu thang điểm

- 1 điểm: ellipse và điểm truy vấn ngẫu nhiên, được vẽ ngay khi chạy.
- 3 điểm: phương pháp Lagrange/đa thức bậc bốn; chạy khi nhấn Enter/Space.
- 1 điểm: vẽ hai điểm cực trị và hai đoạn nối tới điểm truy vấn.

Trong ZIP đã có mã nguồn và ảnh kết quả thực sự xuất từ chương trình.
Đề dùng từ “screenshot”: để bám sát yêu cầu, nên chạy trên máy, nhấn Enter,
chụp màn hình cửa sổ và thay/thêm ảnh đó vào ZIP trước khi nộp.
Đọc và chạy thử để nắm được cách làm, nhất là phần Lagrange và xử lý nghiệm.
