"""Add pipeline table

Revision ID: cf4962b43209
Revises: aba61a88be54
Create Date: 2020-09-30 11:25:37.650804

"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import JSONB


# revision identifiers, used by Alembic.
revision = 'cf4962b43209'
down_revision = 'aba61a88be54'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'pipeline',
        sa.Column('id', sa.Integer(), primary_key=True, index=True),
        sa.Column('project_id', sa.Integer(), sa.ForeignKey("project.id", ondelete="CASCADE"), index=True, nullable=False),
        sa.Column('name', sa.String()),
        sa.Column('description', sa.String()),
        sa.Column('steps', JSONB()),
        sa.Column('created_at', sa.DateTime(), default=sa.sql.func.now(), nullable=False),
    )


def downgrade():
    op.drop_table('pipeline')
