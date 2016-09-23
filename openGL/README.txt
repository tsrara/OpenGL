******인풋 파일 순서******

1.점 찍는 크기
2.색칠할 좌표들 최소값 -> Resolution 변경 가능

3.구의 개수

반지름
중점의 위치
Ka, Kd, Ks, n
반사계수, 굴절계수, 굴절 index(n1/n2)

4.삼각형 개수
좌표 3개
Ka, Kd, Ks, n
반사계수, 굴절계수, 굴절 index(n1/n2)

//gluLookat() 구현을 위한
5. eye position -> Viewpoint 변경 가능
6. look at vector
7. up vector

//Frustum을 위한
8. Near plane
9. Far plane
10. 윈도우 크기(너비,높이)

11. 광원의 개수
position, 빛의 세기

12. attenuation을 위한 상수/일차/이차항 값

13. ambient 빛의 세기


******키보드 버튼******
q,w : eye position x좌표 +/-
a,s : eye position y좌표 +/-
z,x : eye position z좌표 +/-

p : 0번 light position x좌표 +
o : 0번 light position y좌표 +
l : 0번 light position z좌표 +
