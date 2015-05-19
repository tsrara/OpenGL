#define _CRT_SECURE_NO_WARNINGS


#include <OpenGL/gl.h>
#include <GLUT/glut.h>
#include <stdio.h>
#include <vector>
#include <math.h>
#include <list>

using namespace std;


class OctreeNode {
    enum eValue { Depth_Limit = 4, Width = 1000}; //노드의 세부수준을 4단계로 제한, 최초 노드의 가로 크기
    BoundingBox*    m_pBoundBox;                   //현재 노드의 바운딩 박스
    list <Object* >  m_listObjects;                //현재 노드에 포함되는 정적 오브젝트 목록
    bool m_bCulled;                                //이 노드가 컬링 되었는지 여부
    OctreeNode*     m_pParent;                     //부모 노드
};

struct Vector3{
    float x, y, z;
}vec;

void AddChildNode (OctreeNode* );                  //이 노드에 자식 노드를 추가합니다.
bool AddObject (Object* );                         //이 노드에 오브젝트를 추가합니다
void CullNode (bool);                              //이 노드를 컬링합니다.
OctreeNode* const GetChildNode(size_t);            //자식 노드를 가져옵니다.
void SetParent(OctreeNode* const);                 //부모 노드를 설정합니다.
OctreeNode* const GetParent();                     //부모 노드를 반환합니다.
bool IsInNode(const GLTVector3);                   //현재 노드에 해당 정점이 포함되는지 판단합니다.
void Render();                                     //이 노드의 물체들을 렌더링합니다.
bool IsCulled();                                   //이 노드가 컬링되었는지 판단합니다.

//vCenter : 노드의 중심 위치, fHalfWidth : 노드의 반지름 depthLimit : 현재 노드 진입 단계
OctreeNode* BuildOctree(Vector3 vCenter, float fHalfWidth, int depthLimit){
    
    //제한된 진입단계에 도달하면 더이상 자식노드를 생성하지 않습니다.
    if(depthLimit < 0){
        return NULL;
    }//if
    
    //현재 노드를 생성
    OctreeNode* pOctNode = new OctreeNode();
    BoudingBox* pBBox = pOctNode->GetBoundingBox();
    pOctNode->SetPosition(vCenter);
    pBBox->SetRadius(fHalfWidth);
    
    //재귀적으로 8개의 자식노드들을 생성합니다.
    GLTVector3 vOffset;
    GLTVector3 vChildCenter;
    float fStep = fHalfWidth * 0.5f;
    
    //8개의 자식 노드들에 대해서 중심 위치를 설정하고 트리를 생성
    for (int iTree=0; iTree < OctreeNode :: Child_Node_Count; ++iTree){
        vOffset[0] = ((iTree & 1)? fStep : -fStep);
        vOffset[1] = ((iTree & 4)? fStep : -fStep);
        vOffset[2] = ((iTree & 2)? fStep : -fStep);
        vChildCenter[0] = vOffset[0] + vCenter[0];
        vChildCenter[1] = vOffset[1] + vCenter[1];
        vChildCenter[2] = vOffset[2] + vCenter[2];
        pOctNode->AddChildNode( BuildOctree( vChildCenter, fStep, depthLimit - 1 ) );
    }//for
    
    return pOctNode;
}//void BuildOctree(float fCenter, float fHalfWidth, size_t Depth)



OctreeNode* pOctreeNode = BuildOctree( Vector3( 0, 0, 0 ), 1000, 4 );
 //최상위 노드의 중심 위치, 최상위 노드의 반지름, 자식 노드 생성의 제한 레벨 )

OctreeNode* pSecondNode = pOctreeNode->GetChildrenNode( 2 ); //8개의 자식노드 중 2번의 자식 노드를 가져온다.
 
//vObjPos : 판단하고자 하는 물체의 위치,    pNode : 판단하고자 하는 노드
bool FindCurrentPosNode(const Vector3 vObjPos, OctreeNode* const pNode){
    
    Vector3* pvNodePos = pNode->GetPosition();
    
    float fRadius = pNode->GetBoundingBox()->GetRadius();
    //fRadius *= LOOSE_FACTOR //느슨한 옥트리 같은 경우 반지름에 계수를 곱해줄 수 있다.
    float fMin, fMax;
    
    //이 부분은 현재 물체가 이 노드 안에 완전히 포함되는지 판단합니다.
    for(int index = 0; index < 3; ++index){
        fMin = (*pvNodePos)[index] - fRadius;
        fMax = (*pvNodePos)[index] + fRadius;
        
        //포함되지 않으면 실패를 반환하죠
        if(vObjPos[index] < fMin || vObjPos[index] > fMax){
            return false;
        }
        
        //만약 물체가 노드 안에 완전히 속한다면 해당 노드를 반환, 자식 노드들도 판단해 보기
        return true;
}


void display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    glEnable(GL_DEPTH_TEST);
    // view transform
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    gluLookAt(0, 0, 10, 0, 0, 0, 0, 1, 0);
    glutSwapBuffers();
    glutPostRedisplay();
}

void reshape(int w, int h)
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    //gluPerspective(30.0, w/h, 1, 100);
    glOrtho(-1.5, 1.5, -0.3, 2.2, 0, 100);
    glutPostRedisplay();
}

 
 int main(int argc, char** argv)
 {
    glutInit(&argc, argv);
 
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowSize(800, 600);
    glutCreateWindow("2012210067_term");
    glClearColor(0.0, 0.0, 0.0, 1.0);
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    
    glutMainLoop();
    
    return 0;
 }

